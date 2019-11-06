import math
from sympy import *
import sympy
import scipy.optimize
import numpy as np

'''
Hanwei Wangd RHDHV 15_2_2019
'''


w_wave = 0.5
water_depth = 10


def f3(x):

    return  9.81 * (1/0.1)*(1/0.1)/2/math.pi * math.tanh(2*math.pi*10/x) - x


solve_3 = scipy.optimize.fsolve(f3, [1])
print(solve_3)

cor_x = range(12) # remember the length of x is always 2 more than the real case
cor_y = [3]*11

def simpr(x, y, n):
    s1 = 0
    s2 = 0
    h = (x[-1] - x[0])/(2*len(x))
    index_1 = np.arange(x[0], x[-2])
    index_2 = np.arange(x[1], x[-1])

    for index in index_1:
        s1 += y[index]
        print(index)
    for index in index_2:
        s2 = s2+ y[index]
    print('{}, {}'.format(s1, s2))
    s = h * (y[0] + y[-1] + 4*s1 + 2*s2)/3
    return s

print(simpr(cor_x, cor_y, 10))


