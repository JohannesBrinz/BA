import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath


#Defining important constants

N = int(4)
M = int(8)


#Defining matricies
mathcalA = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
kappa = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
U = np.zeros(shape=(N*N,N*N))
tildekappa = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
K_0 = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
kappa_0 = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
mathcalL = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
Lambda = np.array([])
dif = np.array([])


#Implementing U from Timm and Lange
for i in range(N):
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if (N*m+n) == (N*i +j):
                    if m == n:
                        U[N*i+j][(N*m+n)] = 1/np.sqrt((m+2) * (m+1))
                    else:
                        U[N*i+j][(N*m+n)] = 1
                elif n == m:
                    if j == i:
                        if n<i:
                            U[N*i+j][(N*m+n)] = 1/np.sqrt((m+2) * (m+1))
                        elif n == i+1:
                            U[N*i+j][(N*m+n)] = -(m+1)/np.sqrt((m+2) * (m+1))


#Create \mathcalA
for i in range(N):
    for j in range(N):
        for m in range(N):
            for n in range(N):
                mathcalA[(N*i+j)][(N*m+n)] = 13

for i in range(N):                  #Block 1
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if (N*i+j) < int(N*N/4):
                    if int(N*N/4) <= (N*m+n):
                        if (N*m+n) < int(3*N*N/4):
                            mathcalA[(N*i+j)][(N*m+n)] = 0

for i in range(N):                  #Block 2
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if (N*m+n) < int(N*N/4):
                    if int(N*N/4) <= (N*i+j):
                        if (N*i+j) < int(3*N*N/4):
                            mathcalA[(N*i+j)][(N*m+n)] = 0

for i in range(N):                  #Block 3
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if (N*m+n) >= int(3*N*N/4):
                    if int(N*N/4) <= (N*i+j):
                        if (N*i+j) < int(3*N*N/4):
                            mathcalA[(N*i+j)][(N*m+n)] = 0

for i in range(N):                  #Block 4
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if (N*i+j) >= int(3*N*N/4):
                    if int(N*N/4) <= (N*m+n):
                        if (N*m+n) < int(3*N*N/4):
                            mathcalA[(N*i+j)][(N*m+n)] = 0


for i in range(N):                  #A_ijmn = A_jinm
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if mathcalA[(N*i+j)][(N*m+n)] == 0:
                    mathcalA[(N*j+i)][(N*n+m)] = mathcalA[(N*i+j)][(N*m+n)]

plt.matshow(mathcalA.real)
plt.colorbar()
plt.show()


'''for i in range(N):                  #Normal distributed entrys sigma = 1 for diagonal, sigma = 1/2 off diagonal
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if mathcalA[(N*i+j)][(N*m+n)] == 13:
                    if m == n:            #real elements
                        if i==j:
                            mathcalA[(N*i+j)][(N*m+n)] = np.random.normal(loc=0.0, scale=1, size=None)
                            mathcalA[(N*j+i)][(N*n+m)] = mathcalA[(N*i+j)][(N*m+n)]

                        #complex entries
                        else:
                            mathcalA[(N*i+j)][(N*m+n)] = complex(np.random.normal(loc=0.0, scale=0.5, size=None), np.random.normal(loc=0.0, scale=1, size=None))
                            mathcalA[(N*j+i)][(N*n+m)] = np.conj(mathcalA[(N*i+j)][(N*m+n)])
                    else:
                        mathcalA[(N*i+j)][(N*m+n)] = complex(np.random.normal(loc=0.0, scale=0.5, size=None), np.random.normal(loc=0.0, scale=0.5, size=None))
                        mathcalA[(N*j+i)][(N*n+m)] = np.conj(mathcalA[(N*i+j)][(N*m+n)])
'''
#kappa = AA^T
kappa = mathcalA.dot(mathcalA.T)

#tildekappa = UkappaU^T
tildekappa = U.dot(kappa.dot(U.T))

K_0 = tildekappa
for a in range(N*N):
        K_0[N*(N-1)+(N-1)][a] = 0
for b in range(N*N):
        K_0[b][N*(N-1)+(N-1)] = 0


#Calculating mathcalL from K_0 via Lange,Timm (easy way, by back-transforming) Eq.:(23) (bad behavior)
kappa_0 = U.T.dot(K_0.dot(U))

for i in range(N):
    for j in range(N):
        for m in range(N):
            for n in range(N):
                if i == m:
                    for sum in range(N):
                        mathcalL[N*i+j][N*m+n] += 0.5*kappa_0[sum*N+j][sum*N+n]
                if j == n:
                    for sum in range(N):
                        mathcalL[N*i+j][N*m+n] += 0.5*kappa_0[i*N+sum][m*N+sum]
                mathcalL[N*i+j][N*m+n] += kappa_0[N*i+m][N*j+n]

values = lg.eigvals(mathcalL)

Lambda = np.append(Lambda, values)
