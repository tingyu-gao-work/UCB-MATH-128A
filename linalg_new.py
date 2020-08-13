# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 14:05:34 2020

@author: 高廷宇
"""

import numpy as np
class lg_solver():
    @staticmethod
    def __backward(U, b):
        m, n = np.shape(U)
        i = m
        while i > 0:
            i += -1
            s = 0
            for j in range(i+1, m, 1):
                s += b[j] * U[i, j]
            b[i] = 1/U[i, i] * (b[i] - s)
        return b
    
    @staticmethod
    def __forward(L, b):
        m, n = np.shape(L)
        i = m
        for i in range(m):
            s = 0
            for j in range(i):
                s += b[j] * L[i, j]
            b[i] = 1/L[i, i] * (b[i] - s)
        return b   
        
    @staticmethod
    def __change_row(A, s):
        # Default: Partial Pivoting
        interchange = False # Initial value for whether line changes
        
        B = np.copy(A)
        first_col = B[:,0]
        
        if s == 'SPP': # Scaled Partial Pivoting
            m, n = np.shape(B)
            for i in range(m):
                # Find the largest value each line
                largest = np.amax(np.absolute(B[i,:]))
                if largest != 0:
                    B[i,:] = 1/largest * B[i,:]  
            first_col = B[:,0]
        
        # After renormalization, find the largest element in the first col
        index = np.argmax(np.absolute(first_col))
        
        P = [] # Permutation (initial value)
        
        if index != 0:
            # Change rows
            A[[0, index]] = A[[index, 0]]
            interchange = True
            
            # Create coorresponding P
            P = [0, index]
            
        return A, interchange, P
    
    @staticmethod
    def __gausselim(A, s):
        interchange = False # Initial Value for whether row changes
        m, n = np.shape(A)
        P = np.eye(m) #Permutation matrix
        M = np.zeros((m,m)) #Operation
        
        if s is None:
            if A[0,0] == 0:
                A, interchange, P = lg_solver.__change_row(A, s)
        elif s == 'PP' or s == 'SPP':
            A, interchange, P = lg_solver.__change_row(A, s)
        else:
            raise ValueError('No this strategy.')
        
        for i in range(1, m):
            M[i,0] = - A[i,0]/A[0,0]
            A[i,:] = A[i,:] + M[i, 0] * A[0,:]
        
        return A[0,:], A[1:,1:], M, interchange, P
    
    def __return_max(self, A):
        m, n = np.shape(A)
        max_list = np.array([np.amax(np.absolute(A[i,:]))\
                             for i in range(m)])
        if np.amin(max_list) == 0:
            raise TypeError('Scaled Partial Pivoting cannot be used.')
        return max_list

    def LU_decomposition(self, A, s = None):
        A = A.astype(float)
        
        max_list = np.zeros(1)
        if s == 'SPP': # Scaled Partial Pivoting
            B = np.copy(A)
            max_list = self.__return_max(A)
            A = np.diag(max_list) @ B
        
        Permutation = []
        temp_L = np.zeros((np.shape(A)[0],np.shape(A)[0]))
        aug = []
        
        compensation = 0
        while np.shape(A)[0] != 0:
            line, A, M, interchange, P = lg_solver.__gausselim(A, s)
            
            if interchange:
                P = (np.array(P) + compensation).tolist()
                Permutation.append(P)
            
            temp = M[1:, 0]
            row, col = np.shape(M)
            M = np.zeros((row + compensation, col + compensation))
            M[(compensation + 1):, compensation] = temp
            temp_L += M
            
            aug.append(np.hstack((np.zeros(compensation), line)))
            compensation += 1
        
        I = np.eye(np.shape(M)[0])
        L = np.copy(I) - temp_L
        P = np.copy(I)
        n = len(Permutation)%2
        U = np.vstack(aug)
        
        if s == 'SPP':
            U = np.diag(1/max_list) @ U
        
        for index in Permutation:
            P[index] = P[index[::-1]]
            
        return P, n, L, U
    
    def LU_Factorization(self, A): 
        if A[0,0] == 0:
            raise ValueError('Factorization impossible')
        
        m, n = np.shape(A)
        if m!=n:
            raise TypeError('Need Square Matrix')
        
        L = np.eye(m)
        U = np.zeros((m,n))
        U[0,0] = A[0,0]
        
        for j in range(1,n):
            U[0,j] = A[0,j]
            L[j,0] = A[j,0]/U[0,0]
            
        for i in range(1,n-1):
            s = L[i,0] * U[0,i]
            for k in range(1,i):
                s += L[i,k] * U[k,i]
                
            U[i,i] = A[i,i] - s
            
            if U[i,i]==0.0:
                raise ValueError('Factorization impossible')
            
              
            for j in range(i+1, n):
                s = L[i,0] * U[0,j]
                for k in range(1,i):
                    s += L[i,k] * U[k,j]
                    
                U[i,j] = A[i,j] - s
                
                s = L[j,0] * U[0,i]
                for k in range(1,i):
                    s += L[j,k] * U[k,i]
                    
                L[j,i] = 1/U[i,i] * (A[j, i] - s)
        
        s = L[n-1,0] * U[0,n-1]
        for k in range(1,n-1):
            s += L[n-1,k] * U[k,n-1]
                
        U[n-1,n-1] = A[n - 1, n - 1] - s
        
        if U[n-1,n-1]==0.0:
            raise ValueError('Factorization impossible') 
        return L, U
    
    def LDL_Factorization(self, A):
        L, U = self.LU_Factorization(A)
        return L, np.diag(np.diag(U))
        
    
    def det(self, A, s = None):
        m, n = np.shape(A)
        if m != n:
            raise ValueError('The matrix is not a square matrix.')
        P, n, L, U = self.LU_decomposition(A, s)
        if n == 0:
            return np.prod(np.diag(U))
        return -np.prod(np.diag(U))
    
    def solve(self, A, b, s = None):
        P, n, L, U = self.LU_decomposition(A, s)
        if np.amin(np.absolute(np.diag(U))) == 0.0:
            raise ValueError('The matrix is singular.')       
        y = self.__forward(L, P @ b)
        x = self.__backward(U, y)
        return x
    
    def solve_by_LU(self, A, b):
        L, U = self.LU_Factorization(A)
        if np.amin(np.absolute(np.diag(U))) == 0.0:
            raise ValueError('The matrix is singular.')       
        y = self.__forward(L, b.astype(float))
        x = self.__backward(U, y)
        return x