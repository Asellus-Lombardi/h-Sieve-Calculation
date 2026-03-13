# ここにさっきのPythonコードを貼る
import numpy as np
from scipy.linalg import eigh
from math import factorial
def get_primes(limit):
    """ y 以下の素数をリストアップする (エラトステネスの篩) """
    primes = []
    is_prime = [True] * (limit + 1)
    for p in range(2, limit + 1):
        if is_prime[p]:
            primes.append(p)
            for i in range(p * p, limit + 1, p):
                is_prime[i] = False
    return primes
def calculate_boost(y):
    """ ブースト因子 B(y) = prod(2p/(2p-1)) を計算する """
    primes = get_primes(y)
    boost = 1.0
    for p in primes:
        boost *= (2 * p) / (2 * p - 1)
    return boost
def poly_integral_simplex(a, b):
    """ 多項式 t1^a * t2^b のシンプレックス(t1+t2 <= 1)上での積分公式 """
    # Integral = a! b! / (a+b+2)!
    return factorial(a) * factorial(b) / factorial(a + b + 2)
def construct_matrices(D):
    """ 対称多項式基底を用いて行列 A と M を構成する (k=2) """
    # 基底: phi = t1^a * t2^b + t1^b * t2^a (a <= b, a+b <= D)
    basis = []
    for b in range(D + 1):
        for a in range(b + 1):
            if a + b <= D:
                basis.append((a, b))
    
    L = len(basis)
    A = np.zeros((L, L))
    M = np.zeros((L, L))
    
    for i in range(L):
        for j in range(L):
            a1, b1 = basis[i]
            a2, b2 = basis[j]
            
            # 行列 A (I(F)用) の積分計算
            # (t1^a1*t2^b1 + t1^b1*t2^a1) * (t1^a2*t2^b2 + t1^b2*t2^a2)
            # 展開して4つの項を simplex 積分
            A[i, j] = (poly_integral_simplex(a1+a2, b1+b2) + 
                       poly_integral_simplex(a1+b2, b1+a2) + 
                       poly_integral_simplex(b1+a2, a1+b2) + 
                       poly_integral_simplex(b1+b2, a1+a2))
            
            # 行列 M (J_m(F)用) の積分計算 (m=1, 2 の和)
            # J_1 = integral_0^1 (integral_0^{1-t1} phi dt2)^2 dt1
            # 多項式 F に対する Maynard の解析的な J 積分の形を利用
            def j_term(a_i, b_i, a_j, b_j):
                return (1.0 / ((b_i + 1) * (b_j + 1)) * 
                        factorial(a_i + a_j) * factorial(b_i + b_j + 2) / 
                        factorial(a_i + a_j + b_i + b_j + 3))
            
            # 基底 phi の対称性を考慮
            m_val = (j_term(a1, b1, a2, b2) + j_term(a1, b1, b2, a2) +
                     j_term(b1, a1, a2, b2) + j_term(b1, a1, b2, a2))
            M[i, j] = 2 * m_val # m=1 と m=2 は対称なため 2倍
            
    return A, M
def verify_all_thresholds():
    """ 論文用の検証テーブルを出力する """
    print(f"{'Threshold y':<15} | {'Boost B(y)':<12} | {'Sieve Ratio rho_2(h)':<20} | {'rho > 4?'}")
    print("-" * 65)
    
    # D=6 で rho_2(1) を固定算出
    A, M = construct_matrices(6)
    eigvals = eigh(M, A, eigvals_only=True)
    rho2_1 = np.max(eigvals)
    
    thresholds = [100, 241, 263, 1000, 10000, 30000]
    for y in thresholds:
        b_y = calculate_boost(y)
        rho2_h = rho2_1 * b_y
        passed = "Yes" if rho2_h > 4.0 else "No"
        print(f"{y:<15} | {b_y:<12.6f} | {rho2_h:<20.6f} | {passed}")
if __name__ == "__main__":
    # 多項式次数 D による収束性の確認
    print("Maynard ratio rho_2(1) convergence check:")
    for d in [2, 4, 6, 8]:
        A_d, M_d = construct_matrices(d)
        eigvals_d = eigh(M_d, A_d, eigvals_only=True)
        print(f"D = {d}: rho_2(1) = {np.max(eigvals_d):.6f}")
    print("\n" + "="*65)
    
    # 閾値テーブルの生成
    verify_all_thresholds()