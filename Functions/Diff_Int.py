# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 07:53:55 2020

@author: 高廷宇
"""
import numpy as np
class Diff_Int():
    def Richardson(self, N, j, h):  
        if j == 1:
            return N(h)
        else:
            value1 = self.Richardson(N, j-1, h/2)
            value2 = self.Richardson(N, j-1, h)
            return value1 + (value1 - value2)/(2**(j-1) - 1)

    def Trapezoidal(self, f, a, b):
        h = b-a
        return h/2 * (f(a) + f(b))

    def Composite_Trapezoidal(self, f, a, b, n):
        h = (b-a)/n
        j = np.array(range(1,n,1))
        x = a * np.ones(n-1) + h * j
    
        return h/2 * (f(a) + f(b)) + h * np.sum(np.vectorize(f)(x))

    def Simpson(self, f, a, b):
        h = (b-a)/2
        x0 = a
        x1 = a + h
        x2 = b
        return h/3 * (f(x0) + 4 * f(x1) + f(x2))

    def Composite_Simpson(self, f, a, b, n):
        if n%2 != 0:
            raise ValueError()
    
        h = (b-a)/n
        j = np.array(range(n+1))
        x = a * np.ones(n+1) + h * j
    
        F = np.vectorize(f)
        return h/3 * (f(a) + 2 * np.sum(F(x[2:n:2]))\
                      + 4 * np.sum(F(x[1:n:2])) + f(b))
    
    @staticmethod
    def __Legendre(n):
        from Interpolation import comb
    
        c = np.zeros(n+1)
        for k in range(n+1):
            c[k] = comb(n,k) * comb((n+k-1)/2, n)
        return 2**n * c[::-1]

    def gaussquad(self, n):
        from itertools import combinations as C
    
        x = np.sort(np.roots(self.__Legendre(n)))
        coi = np.zeros((n,n))
    
        for i in range(n):
            for k in range(n):
                l = list(range(n))
                l.remove(i)
            
                combination = C(l, n-1-k)
            
                s = 0
                for tup in combination:
                    indices = np.array(tup)
                    if indices.size == 0:
                        s += 1
                    else:
                        s += np.prod(-x[indices])
                    
                coi[i,k] = s/np.prod(x[i] * np.ones(n-1) - x[np.array(l)])
    
        def temp(A):
            s = 0
            for i in range(0,len(A),2):
                s += 2/(i+1) * A[i]
            return s
    
        c = np.zeros(n)
        for i in range(n):
            c[i] = temp(coi[i,:])
        return x, c

    def Gaussian_Quadrature(self, f, a, b, n):
        x, c = self.gaussquad(n)
    
        def F(x):
            return (b-a)/2 * f(((b-a)*x + b+a)/2)
    
        return np.sum(c * np.vectorize(F)(x))
