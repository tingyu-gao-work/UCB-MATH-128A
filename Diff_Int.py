# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 07:53:55 2020

@author: 高廷宇
"""


def Richardson(N, j, h):  
    if j == 1:
        return N(h)
    else:
        value1 = Richardson(N, j-1, h/2)
        value2 = Richardson(N, j-1, h)
        return value1 + (value1 - value2)/(2**(j-1) - 1)

def Trapezoidal(f, a, b):
    h = b-a
    return h/2 * (f(a) + f(b))

def Simpson(f, a, b):
    h = (b-a)/2
    x0 = a
    x1 = a + h
    x2 = b
    return h/3 * (f(x0) + 4 * f(x1) + f(x2))