#!/home/daiki/keio/.venv/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 
# データ列から7次関数の係数を抽出
data = "2,Mean,235.771,4.72732e-10,-3.05675,4.72733e-10,0.247008,4.72745e-10,0.0178393,4.73129e-10,0.000937921,4.78614e-10,-0.000146824,1.44331e-09,4.13281e-05,0.0720827,5.18158e-08,1.36741e-11,591746,10,228.895,3.93189"
values = data.split(',')
# 7次関数の係数は、3番目から17番目の
# 奇数番目の値に対応
coefficients = [float(values[i]) for i in range(2, 18, 2)]
# 7次関数を定義
def polynomial_7th(x, coeffs):
    return sum(c * x**i for i, c in enumerate(coeffs))
# xの範囲を設定
x = np.linspace(0, 900, 400)
# 7次関数のy値を計算
y = polynomial_7th(x, coefficients)
# プロット
plt.figure(figsize=(10, 6))
plt.plot(x, y, label='7th Degree Polynomial Fit', color='blue')
plt.title('7th Degree Polynomial Fit from Data')
plt.xlabel('x')
plt.ylabel('f(x)')                  
plt.legend()
plt.grid()
plt.show()