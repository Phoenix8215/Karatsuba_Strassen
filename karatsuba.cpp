#include <bits/stdc++.h>
using namespace std;

/*
 * 将字符串转换为逆序存储的 vector<int>，如 "12345" -> [5,4,3,2,1]
 * 方便低位对齐进行大数运算
 */
vector<int> stringToVector(const string &num) {
  vector<int> result(num.size());
  for (int i = 0; i < (int)num.size(); i++) {
    result[num.size() - 1 - i] = num[i] - '0';
  }
  return result;
}

/*
 * 去除字符串表示的数字前导零（最高位的 0）
 */
string removeLeadingZeros(const string &str) {
  int pos = 0;
  while (pos < (int)str.size() - 1 && str[pos] == '0') {
    pos++;
  }
  return str.substr(pos);
}

/*
 * 将逆序 vector<int> 转回字符串（去掉前导零，若全0则保留一个0）
 */
string vectorToString(const vector<int> &v) {
  string result;
  for (int i = (int)v.size() - 1; i >= 0; i--) {
    result.push_back((char)(v[i] + '0'));
  }
  return removeLeadingZeros(result);
}

/*
 * 常规大整数乘法 (竖式)
 * 假设 A 和 B 都是逆序存储的数字，返回值也为逆序存储
 */
vector<int> NaiveMultiply(const vector<int> &A, const vector<int> &B) {
  vector<int> product(A.size() + B.size(), 0);
  for (int i = 0; i < (int)A.size(); i++) {
    long long carry = 0;
    for (int j = 0; j < (int)B.size() || carry; j++) {
      long long sum =
          product[i + j] + 1LL * A[i] * (j < (int)B.size() ? B[j] : 0) + carry;
      product[i + j] = (int)(sum % 10);
      carry = sum / 10;
    }
  }
  while (product.size() > 1 && product.back() == 0)
    product.pop_back();
  return product;
}

/*
 * 对长度不足 n 的向量进行补齐（低位对齐）
 * 方便 Karatsuba 递归处理时拆分
 */
void padZeros(vector<int> &v, int n) {
  while ((int)v.size() < n) {
    v.push_back(0);
  }
}

/*
 * 对大数进行加法：C = A + B
 * A, B 均为逆序存储
 */
vector<int> addVec(const vector<int> &A, const vector<int> &B) {
  int n = max(A.size(), B.size());
  vector<int> result;
  result.reserve(n + 1);

  long long carry = 0;
  for (int i = 0; i < n || carry; i++) {
    long long x = (i < (int)A.size() ? A[i] : 0);
    long long y = (i < (int)B.size() ? B[i] : 0);
    long long sum = x + y + carry;
    result.push_back(sum % 10);
    carry = sum / 10;
  }
  return result;
}

/*
 * 对大数进行减法：C = A - B
 * 假设 A >= B，A, B 均为逆序存储，返回值也为逆序存储
 */
vector<int> subVec(const vector<int> &A, const vector<int> &B) {
  // A >= B
  vector<int> result;
  result.reserve(A.size());

  long long carry = 0;
  for (int i = 0; i < (int)A.size(); i++) {
    long long x = A[i] - carry;
    long long y = (i < (int)B.size() ? B[i] : 0);
    if (x < y) {
      x += 10;
      carry = 1;
    } else {
      carry = 0;
    }
    result.push_back((int)(x - y));
  }
  while (result.size() > 1 && result.back() == 0) {
    result.pop_back();
  }
  return result;
}

/*
 * Karatsuba 乘法核心函数，对逆序存储的大数 A, B 进行乘法
 * 返回结果同样是逆序存储
 */
vector<int> Karatsuba(const vector<int> &A, const vector<int> &B) {
  int n = max(A.size(), B.size());
  // 基础情形：如果足够小，直接使用常规乘法
  if (n <= 64) {
    return NaiveMultiply(A, B);
  }

  // 长度补齐
  vector<int> A_(A), B_(B);
  padZeros(A_, n);
  padZeros(B_, n);

  int half = n / 2;

  // 拆分为 A = A1 * 10^(half) + A0，B = B1 * 10^(half) + B0
  vector<int> A0(A_.begin(), A_.begin() + half);
  vector<int> A1(A_.begin() + half, A_.end());
  vector<int> B0(B_.begin(), B_.begin() + half);
  vector<int> B1(B_.begin() + half, B_.end());

  // 递归计算
  // z2 = A1 * B1
  vector<int> z2 = Karatsuba(A1, B1);

  // z0 = A0 * B0
  vector<int> z0 = Karatsuba(A0, B0);

  // z1 = (A1 + A0)*(B1 + B0) - z2 - z0
  vector<int> A1A0 = addVec(A1, A0);
  vector<int> B1B0 = addVec(B1, B0);
  vector<int> z1 = Karatsuba(A1A0, B1B0);
  z1 = subVec(z1, z2);
  z1 = subVec(z1, z0);

  // 将结果合并
  // result = z2 * 10^(2*half) + z1 * 10^half + z0
  // 对应逆序存储，需要在合并时左移
  // 左移 half（乘以 10^half）
  auto shiftLeft = [&](vector<int> &v, int k) { v.insert(v.begin(), k, 0); };

  shiftLeft(z2, 2 * half);
  shiftLeft(z1, half);

  vector<int> result = addVec(z2, z1);
  result = addVec(result, z0);

  // 去掉末端多余的 0
  while (result.size() > 1 && result.back() == 0) {
    result.pop_back();
  }

  return result;
}

/*
 * 对外暴露的 Karatsuba 大整数乘法接口
 * 输入和输出均为 string
 */
string KaratsubaMultiply(const string &num1, const string &num2) {
  if (num1 == "0" || num2 == "0")
    return "0";

  // 转化为逆序存储
  vector<int> A = stringToVector(num1);
  vector<int> B = stringToVector(num2);

  vector<int> result = Karatsuba(A, B);

  return vectorToString(result);
}

/*
 * 对外暴露的 Naive 大整数乘法接口
 * 输入和输出均为 string
 */
string NaiveMultiply(const string &num1, const string &num2) {
  if (num1 == "0" || num2 == "0")
    return "0";

  vector<int> A = stringToVector(num1);
  vector<int> B = stringToVector(num2);

  vector<int> product = NaiveMultiply(A, B);

  return vectorToString(product);
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  const int TEST_TIMES = 5;

  // 测试数字的长度（位数）
  // 可以修改以观察在不同规模下的性能差异
  vector<int> lengthList = {1000, 2000, 4000, 8000, 16000};

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> dist(0, 9);

  for (auto len : lengthList) {
    cout << "数字位数: " << len << "\n";

    long long totalNaiveTime = 0;
    long long totalKaratsubaTime = 0;

    for (int t = 0; t < TEST_TIMES; t++) {
      // 生成两个随机大整数
      string num1, num2;
      // 确保最高位不为 0
      num1.push_back('1' + dist(gen) % 9);
      num2.push_back('1' + dist(gen) % 9);

      for (int i = 1; i < len; i++) {
        num1.push_back((char)('0' + dist(gen)));
        num2.push_back((char)('0' + dist(gen)));
      }

      // 记录 NaiveMultiply 耗时
      auto startNaive = chrono::high_resolution_clock::now();
      volatile auto naiveRes = NaiveMultiply(num1, num2);
      auto endNaive = chrono::high_resolution_clock::now();
      long long naiveDuration =
          chrono::duration_cast<chrono::microseconds>(endNaive - startNaive)
              .count();
      totalNaiveTime += naiveDuration;

      // 记录 KaratsubaMultiply 耗时
      auto startKara = chrono::high_resolution_clock::now();
      volatile auto karaRes = KaratsubaMultiply(num1, num2);
      auto endKara = chrono::high_resolution_clock::now();
      long long karaDuration =
          chrono::duration_cast<chrono::microseconds>(endKara - startKara)
              .count();
      totalKaratsubaTime += karaDuration;
    }

    cout << "  Naive 平均耗时(us)       : " << (totalNaiveTime / TEST_TIMES)
         << "\n";
    cout << "  Karatsuba 平均耗时(us)   : " << (totalKaratsubaTime / TEST_TIMES)
         << "\n\n";
  }

  return 0;
}
