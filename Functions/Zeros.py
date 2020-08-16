# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 14:53:08 2020

@author: 高廷宇
"""

def Newton(f, df, p0, tol, table = False):
    if table:
        # Print header
        print(' n         p          |p-p0|   ')
        print('-------------------------------')
        
    p = p0
    n = 0
    while True:
        n += 1
        p_old = p
        p = p - f(p)/df(p)
        if table:
            print('{0:2d}  {1:12.8f}  {2:12.8f}'.format(n, p, abs(p-p_old)))
        if abs(p_old - p) < tol:
            break
    return p

def Modified_Newton(f, df, ddf, p0, tol, table = False):
    # Solve f(p) = 0 using Newton's method.
    
    if table:
        # Print header
        print(' n         p          |p-p0|   ')
        print('-------------------------------')
        
    p = p0
    n = 0
    while True:
        n += 1
        p_old = p
        p = p - (f(p) * df(p))/(df(p)**2 - f(p)*ddf(p))
        if table:
            print('{0:2d}  {1:12.8f}  {2:12.8f}'.format(n, p, abs(p-p_old)))
        if abs(p_old - p) < tol:
            break
    return p

def Secant(f, p0, p1, tol, table = False):
    if table:
        # Print header
        print(' n        p_n         p_{n-1}    |p_n-p_{n-1}|  ')
        print('-----------------------------------------------')
    
    n = 0
    while True:
        n += 1
        sec_line_slope = (f(p0) - f(p1)) / (p0 - p1)
        p_new = p1 - f(p1)/sec_line_slope
        p0 = p1
        p1 = p_new
        if table:
            print('{0:2d}  {1:12.8f}  {2:12.8f}  {3:12.8f}'.format(n, p1, p0, abs(p1-p0)))
        if abs(p1 - p0) < tol:
            break
    return p1
