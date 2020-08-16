# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 20:25:47 2020

@author: 高廷宇
"""


def L(X, k):
    n = len(X) - 1
    def func(x):
        y = 1
        for i in range(n+1):
            if i != k:
                y = y * (x-X[i]) / (X[k] - X[i])
        return y
    return func

def Lagrange(X, f):
    n = len(X) - 1
    
    if len(X) != len(f):
        raise ValueError()
    
    def poly(x):
        y = 0
        for i in range(n+1):
            y += L(X, i)(x) * f[i]
        return y
    return poly



#Combination#
def comb(x, n):
    ans = 1
    for i in range(n):
        ans = ans * ((x - i)/(i+1))
    return ans
#############

#Combination#
def Delta(f, k):
    if k == 0:
        return f[0]
    else:
        f1 = f[0:-1]
        f2 = f[1:]
        g = f2 - f1
        return Delta(g, k-1)
#############

#Combination#
def Nabla(f, k):
    if k == 0:
        return f[-1]
    else:
        f1 = f[0:-1]
        f2 = f[1:]
        g = f2 - f1
        return Nabla(g, k-1)
#############

def Newton_Forward(X0, h, f):
    n = len(f) - 1
    def poly(x):
        s = (x - X0)/h
        y = f[0]
        for k in range(1, n+1, 1):
            temp = comb(s, k) * Delta(f, k)
            y += temp
        return y
    return poly

def Newton_Backward(Xn, h, f):
    n = len(f) - 1
    def poly(x):
        s = (x - Xn)/h
        y = f[-1]
        for k in range(1, n+1, 1):
            if k%2 == 0:
                temp = comb(-s, k) * Nabla(f, k)
                y += temp
            else:
                temp = -comb(-s, k) * Nabla(f, k)
                y += temp
        return y
    return poly


def H(X, j):
    n = len(X) - 1
    
    sum_of_dL = 0
    for i in range(n+1):
        if i!=j:
            sum_of_dL += 1/(X[j] - X[i])
    
    def func(x):
        return (1 - 2 * (x - X[j]) * sum_of_dL) * L(X, j)(x)**2
    return func

def hH(X, j):
    def func(x):
        return (x - X[j]) * L(X, j)(x)**2
    return func

def Hermite(X, f, df):
    n = len(X) - 1
    
    if (len(X) != len(f)) or (len(X) != len(df)) or (len(df) != len(f)):
        raise ValueError()
        
    def poly(x):
        y = 0
        for i in range(n+1):
            y += H(X, i)(x) * f[i] + df[i] * hH(X, i)(x)
        return y
    return poly







def Cubic_Spline(X, f, df = None, kind = 'natural', Return = 'poly'):
    import numpy as np
    
    h = np.array(X[1:]) - np.array(X[0:-1])
    h_0_2 = h[0:-1]
    h_1_1 = h[1:]
    
    a = np.array(f)
    b = np.array([0])
    c = np.copy(b)
    d = np.copy(c)
    
    if kind == 'natural':
        diag = np.insert(np.append((h_0_2 + h_1_1)*2, 1),0,1)
        A = np.diag(np.append(h_0_2, 0), -1)\
            + np.diag(diag, 0)\
                + np.diag(np.insert(h_1_1,0,0), 1)
        b0 = 3 * (a[2:] - a[1:-1]) / h[1:] - 3 * (a[1:-1] - a[0:-2]) / h[0:-1]
        b = np.insert(np.append(b0, 0),0,0)
        
    elif kind == 'clamped':
        if df is None or len(df) != 2:
            raise ValueError()
            
        diag = np.insert(np.append((h_0_2 + h_1_1)*2, 2*h[-1]),0, 2*h[0])
        A = np.diag(h, -1) + np.diag(diag, 0) + np.diag(h, 1)
        
        b0 = 3 * (a[1:] - a[0:-1]) / h
        b1 = np.append(b0, 3 * df[1]) 
        b2 = np.insert(b0, 0, 3*df[0])
        b = b1 - b2
        
    else:
        raise ValueError()
        
    c = np.linalg.solve(A, b)
    b = (a[1:] - a[0:-1])/h - h * (2 * c[0:-1] + c[1:])/3
    d = (-c[0:-1] + c[1:])/(3 * h)
    
    if Return == 'coefficient':
        return a[:-1], b, c[:-1], d
    elif not (Return == 'poly'):
        raise ValueError()
        
    def s(x):
        j = 0
        indi = np.where(X[0:-1] < x)[0]
        if indi.size != 0:
            j = indi[-1]
        return a[j] + b[j]*(x - X[j]) + c[j]*(x - X[j])**2 + d[j]*(x - X[j])**3
    
    return s
