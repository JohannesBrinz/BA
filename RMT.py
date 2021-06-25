import numpy as np
from scipy import linalg as lg
import matplotlib.pyplot as plt
import pandas as pd


#Defining important constants

N = int(4)
M = int(100)


#Defining matricies
S = np.zeros(shape=(N,N))
mathcalS = np.zeros(shape=(N*N,N*N))
mathcalA = np.zeros(shape=(N*N,N*N))
kappa = np.zeros(shape=(N*N,N*N))
U = np.zeros(shape=(N*N,N*N))
tildekappa = np.zeros(shape=(N*N,N*N))
K_0 = np.zeros(shape=(N*N,N*N))
Lambda = []


#Defining S and \mathcalS
for i in range(int(N/2)):
    S[i][i] = 1
for i in range(int(N/2), N):
    S[i][i] = -1

mathcalS = np.kron(S, S)

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

#print(U)

#Generating ensemble

for e in range(M):
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

    for i in range(N):                  #Normal distributed entrys sigma = 1 for diagonal, sigma = 1/2 off diagonal
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    if mathcalA[(N*i+j)][(N*m+n)] == 13:
                        if (N*i+j)==(N*m+n):
                            mathcalA[(N*i+j)][(N*m+n)] = np.random.normal(loc=0.0, scale=1.0, size=None)
                        else:
                            mathcalA[(N*i+j)][(N*m+n)] = np.random.normal(loc=0.0, scale=0.5, size=None)

    #kappa = AA^T
    kappa = mathcalA.dot(mathcalA.T)

    #tildekappa = UkappaU^T
    tildekappa = U.dot(kappa.dot(U.T))

    K_0 = tildekappa
    for a in range(N*N):
            K_0[N*(N-1)+(N-1)][a] = 0
    for b in range(N*N):
            K_0[b][N*(N-1)+(N-1)] = 0

    values, vectors = lg.eig(K_0)

    Lambda.append(values)

print(Lambda)

#plotting histogram
