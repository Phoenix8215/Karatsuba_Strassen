#include <bits/stdc++.h>
using namespace std;

/*
 * 常规（Naive）矩阵乘法：
 * C = A * B, 其中 A, B, C 均为 n x n 矩阵
 * 时间复杂度 O(n^3)
 */
vector<vector<int>> naiveMultiply(const vector<vector<int>> &A,
                                  const vector<vector<int>> &B) {
  int n = A.size();
  vector<vector<int>> C(n, vector<int>(n, 0));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      long long sum = 0;
      for (int k = 0; k < n; k++) {
        sum += 1LL * A[i][k] * B[k][j];
      }
      C[i][j] = (int)sum;
    }
  }
  return C;
}

/*
 * 矩阵加法：C = A + B
 * 要求 A, B, C 同维度
 */
vector<vector<int>> addMatrix(const vector<vector<int>> &A,
                              const vector<vector<int>> &B) {
  int n = A.size();
  vector<vector<int>> C(n, vector<int>(n, 0));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
  return C;
}

/*
 * 矩阵减法：C = A - B
 * 要求 A, B, C 同维度
 */
vector<vector<int>> subMatrix(const vector<vector<int>> &A,
                              const vector<vector<int>> &B) {
  int n = A.size();
  vector<vector<int>> C(n, vector<int>(n, 0));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = A[i][j] - B[i][j];
    }
  }
  return C;
}

/*
 * 将矩阵拆分成四块：
 *   A11, A12
 *   A21, A22
 * 用于 Strassen 递归。
 *
 * in：原矩阵 A，大小为 n x n
 * out：将 A 分成 4 个子矩阵，并返回（每个子矩阵是 n/2 x n/2）
 */
void splitMatrix(const vector<vector<int>> &A, vector<vector<int>> &A11,
                 vector<vector<int>> &A12, vector<vector<int>> &A21,
                 vector<vector<int>> &A22) {
  int n = A.size();
  int mid = n / 2;
  // 分配子矩阵空间
  A11.assign(mid, vector<int>(mid, 0));
  A12.assign(mid, vector<int>(mid, 0));
  A21.assign(mid, vector<int>(mid, 0));
  A22.assign(mid, vector<int>(mid, 0));

  for (int i = 0; i < mid; i++) {
    for (int j = 0; j < mid; j++) {
      A11[i][j] = A[i][j];
      A12[i][j] = A[i][j + mid];
      A21[i][j] = A[i + mid][j];
      A22[i][j] = A[i + mid][j + mid];
    }
  }
}

/*
 * 将 4 块子矩阵合并为一个 n x n 矩阵
 */
void mergeMatrix(vector<vector<int>> &A, const vector<vector<int>> &A11,
                 const vector<vector<int>> &A12, const vector<vector<int>> &A21,
                 const vector<vector<int>> &A22) {
  int mid = A11.size();
  int n = mid * 2;
  A.assign(n, vector<int>(n, 0));
  for (int i = 0; i < mid; i++) {
    for (int j = 0; j < mid; j++) {
      A[i][j] = A11[i][j];
      A[i][j + mid] = A12[i][j];
      A[i + mid][j] = A21[i][j];
      A[i + mid][j + mid] = A22[i][j];
    }
  }
}

/*
 * Strassen 递归乘法核心:
 * 输入: A, B (n x n)
 * 输出: C = A * B (n x n)
 * 要求 n 是 2^k，如果不是，可以补零到最近的 2^k
 */
vector<vector<int>> strassenRecursive(const vector<vector<int>> &A,
                                      const vector<vector<int>> &B) {
  int n = A.size();

  // 对于较小矩阵，直接用 Naive 可减少递归开销
  if (n <= 64) {
    return naiveMultiply(A, B);
  }

  // 拆分为 4 块
  int mid = n / 2;
  vector<vector<int>> A11, A12, A21, A22;
  vector<vector<int>> B11, B12, B21, B22;
  splitMatrix(A, A11, A12, A21, A22);
  splitMatrix(B, B11, B12, B21, B22);

  // M1 = (A11 + A22) * (B11 + B22)
  auto A11A22 = addMatrix(A11, A22);
  auto B11B22 = addMatrix(B11, B22);
  auto M1 = strassenRecursive(A11A22, B11B22);

  // M2 = (A21 + A22) * B11
  auto A21A22 = addMatrix(A21, A22);
  auto M2 = strassenRecursive(A21A22, B11);

  // M3 = A11 * (B12 - B22)
  auto B12B22 = subMatrix(B12, B22);
  auto M3 = strassenRecursive(A11, B12B22);

  // M4 = A22 * (B21 - B11)
  auto B21B11 = subMatrix(B21, B11);
  auto M4 = strassenRecursive(A22, B21B11);

  // M5 = (A11 + A12) * B22
  auto A11A12 = addMatrix(A11, A12);
  auto M5 = strassenRecursive(A11A12, B22);

  // M6 = (A21 - A11) * (B11 + B12)
  auto A21A11 = subMatrix(A21, A11);
  auto B11B12 = addMatrix(B11, B12);
  auto M6 = strassenRecursive(A21A11, B11B12);

  // M7 = (A12 - A22) * (B21 + B22)
  auto A12A22 = subMatrix(A12, A22);
  auto B21B22_ = addMatrix(B21, B22);
  auto M7 = strassenRecursive(A12A22, B21B22_);

  // 组合结果
  // C11 = M1 + M4 - M5 + M7
  auto M1plusM4 = addMatrix(M1, M4);
  auto M1plusM4minusM5 = subMatrix(M1plusM4, M5);
  auto C11 = addMatrix(M1plusM4minusM5, M7);

  // C12 = M3 + M5
  auto C12 = addMatrix(M3, M5);

  // C21 = M2 + M4
  auto C21 = addMatrix(M2, M4);

  // C22 = M1 - M2 + M3 + M6
  auto M1minusM2 = subMatrix(M1, M2);
  auto M1minusM2plusM3 = addMatrix(M1minusM2, M3);
  auto C22 = addMatrix(M1minusM2plusM3, M6);

  // 合并四块
  vector<vector<int>> C;
  mergeMatrix(C, C11, C12, C21, C22);

  return C;
}

/*
 * 若 n 不是 2 的幂，Strassen 通常会将矩阵补零至 2^k x 2^k，再做递归运算。
 * 最后再把结果截回原本尺寸即可。
 */
vector<vector<int>> strassenMultiply(const vector<vector<int>> &A,
                                     const vector<vector<int>> &B) {
  int n = A.size();
  // 假设 A, B 已经是 n x n，且 B.size() == n
  // 1. 找到 >= n 的最小 2^k
  int m = 1;
  while (m < n)
    m <<= 1;

  // 2. 若 m == n，则直接递归，否则补到 m
  if (m == n) {
    return strassenRecursive(A, B);
  } else {
    // 补零
    vector<vector<int>> A_(m, vector<int>(m, 0));
    vector<vector<int>> B_(m, vector<int>(m, 0));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A_[i][j] = A[i][j];
        B_[i][j] = B[i][j];
      }
    }
    // 递归
    auto C_ = strassenRecursive(A_, B_);
    // 截取回 n x n
    vector<vector<int>> C(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        C[i][j] = C_[i][j];
      }
    }
    return C;
  }
}

/*
 * 生成 n x n 矩阵随机数
 */
vector<vector<int>> randomMatrix(int n, int minVal = 0, int maxVal = 10) {
  static random_device rd;
  static mt19937 gen(rd());
  uniform_int_distribution<int> dist(minVal, maxVal);

  vector<vector<int>> mat(n, vector<int>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      mat[i][j] = dist(gen);
    }
  }
  return mat;
}

/*
 * 验证两矩阵是否相等
 */
bool checkEqual(const vector<vector<int>> &X, const vector<vector<int>> &Y) {
  if (X.size() != Y.size())
    return false;
  if (X.size() == 0)
    return true;
  int n = X.size();
  if (X[0].size() != Y[0].size())
    return false;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (X[i][j] != Y[i][j])
        return false;
    }
  }
  return true;
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  // 测试矩阵大小，可根据需要修改(如 128, 256, 512, 1024...)
  vector<int> testSizes = {128, 256, 512, 1024};

  const int TEST_TIMES = 3;

  cout << fixed << setprecision(6);
  for (auto n : testSizes) {
    cout << "========== 矩阵大小 n = " << n << " ==========\n";
    long long totalNaiveTime = 0, totalStrassenTime = 0;

    for (int t = 0; t < TEST_TIMES; t++) {
      // 1. 生成随机矩阵 A, B
      auto A = randomMatrix(n, 0, 9);
      auto B = randomMatrix(n, 0, 9);

      // 2. 测试 NaiveMultiply
      auto startNaive = chrono::high_resolution_clock::now();
      // 移除 volatile 修饰，避免类型不匹配
      auto C_naive = naiveMultiply(A, B);
      auto endNaive = chrono::high_resolution_clock::now();
      long long naiveTime =
          chrono::duration_cast<chrono::microseconds>(endNaive - startNaive)
              .count();
      totalNaiveTime += naiveTime;

      // 3. 测试 StrassenMultiply
      auto startStrassen = chrono::high_resolution_clock::now();
      // 同理，移除 volatile 修饰
      auto C_strassen = strassenMultiply(A, B);
      auto endStrassen = chrono::high_resolution_clock::now();
      long long strassenTime = chrono::duration_cast<chrono::microseconds>(
                                   endStrassen - startStrassen)
                                   .count();
      totalStrassenTime += strassenTime;

      // 4. 检验结果正确性
      if (!checkEqual(C_naive, C_strassen)) {
        cerr << "结果不匹配，可能有错误。\n";
      }
    }

    cout << "Naive 平均耗时(us)     = " << (totalNaiveTime / TEST_TIMES)
         << "\n";
    cout << "Strassen 平均耗时(us) = " << (totalStrassenTime / TEST_TIMES)
         << "\n\n";
  }

  return 0;
}
