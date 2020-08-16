# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 20:56:40 2020

@author: 高廷宇
"""

import numpy as np
class ODE_tools():
    def Euler(self, f, a, b, y0, N):
        # Solve the IVP by using Euler's method
        y0 = np.array([y0]).flatten()
        h = (b - a)/N
        result = np.zeros((N+1, y0.size))
        result[0,:] = y0
        
        for i in range(N):
            result[i+1,:] = result[i,:] + h * f(a + i * h, result[i,:])
        return np.column_stack([np.linspace(a, b, N+1), result])
    
    def rk4(self, f, a, b, y0, N):
        # Solve the IVP by using Runge-Kutta(4) method
        y0 = np.array([y0]).flatten()
        h = (b - a)/N
        result = np.zeros((N+1, y0.size))
        result[0,:] = y0
        
        for i in range(N):
            k1 = h * f(a + i * h, result[i,:])
            k2 = h * f(a + i * h + h/2, result[i,:] + k1/2)
            k3 = h * f(a + i * h + h/2, result[i,:] + k2/2)
            k4 = h * f(a + i * h + h, result[i,:] + k3)
            
            result[i+1,:] = result[i,:] + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return np.column_stack([np.linspace(a, b, N+1), result])
