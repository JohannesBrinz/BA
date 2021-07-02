import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath


#Defining important constants

N = int(2)
M = int(8e3)


#Defining matricies
mathcalA = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
kappa = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
U = np.zeros((N**2, N**2))
tildekappa = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
K_0 = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
kappa_0 = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
mathcalL = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
Lambda = pd.DataFrame([])
dif = pd.DataFrame([])
real_lambda = pd.DataFrame([])
complex_lambda = pd.DataFrame([])
zero_lambda = pd.DataFrame([])


#Implementing U from Timm and Lange
print("\nImplementing U...")
for m in range(N):
    for n in range(N):
        for o in range(N):
            for p in range(N):
                if (N* o + p) == (N*m + n):
                    if n == m:
                        if m == (N-1):
                            U[N*m+n][N*o+p] = 1/np.sqrt(N)
                        else:
                            U[N*m+n][N*o+p] = 1/np.sqrt((m+2) * (m+1))
                    else:
                        U[N*m+n][(N*o+p)] = 1
                elif p == o:
                    if n == m:
                        if m < (N-1):
                            if p<m:
                                U[N*m+n][(N*o+p)] = 1/np.sqrt((m+2) * (m+1))
                            elif p == m+1:
                                U[N*m+n][(N*o+p)] = -(m+1)/np.sqrt((m+2) * (m+1))
                        else:
                            U[N*m+n][(N*o+p)] = 1/np.sqrt(N)


#Generating ensemble
print("\ngenerating ensemble...")
for e in range(M):
    #Create \mathcalA
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    mathcalA[(N*i+j)][(N*m+n)] = complex(np.random.normal(loc=0.0, scale=1, size=None), np.random.normal(loc=0.0, scale=1, size=None)) #=13

    '''for i in range(N):                  #Block 1
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
                        if m == n:            #real elements
                            if i==j:
                                mathcalA[(N*i+j)][(N*m+n)] = np.random.normal(loc=0.0, scale=1, size=None)
                                mathcalA[(N*j+i)][(N*n+m)] = mathcalA[(N*i+j)][(N*m+n)]

                                #complex entries
                            else:
                                mathcalA[(N*i+j)][(N*m+n)] = complex(np.random.normal(loc=0.0, scale=0.5, size=None), np.random.normal(loc=0.0, scale=0.5, size=None))
                                mathcalA[(N*j+i)][(N*n+m)] = np.conj(mathcalA[(N*i+j)][(N*m+n)])
                        else:
                            mathcalA[(N*i+j)][(N*m+n)] = complex(np.random.normal(loc=0.0, scale=0.5, size=None), np.random.normal(loc=0.0, scale=0.5, size=None))
                            mathcalA[(N*j+i)][(N*n+m)] = np.conj(mathcalA[(N*i+j)][(N*m+n)])'''

    #kappa = AA^\dagger
    kappa = mathcalA @ (np.conj(mathcalA).T)

    #tildekappa = UkappaU^T
    tildekappa = np.linalg.multi_dot([U, kappa, U.T])

    K_0 = tildekappa
    for a in range(N*N):
            K_0[N*(N-1)+(N-1)][a] = 0
    for b in range(N*N):
            K_0[b][N*(N-1)+(N-1)] = 0

    #Calculating mathcalL from K_0 via Lange,Timm (easy way, by back-transforming) Eq.:(23)
    kappa_0 = U.T @ K_0 @ U
    mathcalL = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
    for m in range(N):
        for n in range(N):
            for p in range(N):
                for q in range(N):
                    mathcalL[N*m+n][N*p+q] += kappa_0[N*m+p][N*n+q]
                    if p == m:
                        for r in range(N):
                            mathcalL[N*m+n][N*p+q] = mathcalL[N*m+n][N*p+q]-0.5*kappa_0[r*N+n][r*N+q]
                    if q == n:
                        for r in range(N):
                            mathcalL[N*m+n][N*p+q] = mathcalL[N*m+n][N*p+q]-0.5*kappa_0[r*N+p][r*N+m]


    values = lg.eigvals(mathcalL)
    print(values)
    Lambda = np.append(Lambda, values)


#separating real and complex eigenvalues
print("\nseparating complex and real eigenvalues...")
for i in range (len(Lambda)):
    if np.absolute(Lambda[i].real) <= 3e-9:
        if np.absolute(Lambda[i].imag) <= 3e-9:
            zero_lambda = np.append(zero_lambda, Lambda[i])
        else:
            complex_lambda = np.append(complex_lambda, Lambda[i])
    if np.absolute(Lambda[i].real) > 3e-9:
        if np.absolute(Lambda[i].imag) > 3e-9:
            complex_lambda = np.append(complex_lambda, Lambda[i])
        else:
            real_lambda = np.append(real_lambda, Lambda[i])


print("\nNumber of zero eigenvalues: ", len(zero_lambda))
print("\nNumber of real eigenvalues: ", len(real_lambda))
print("\nNumber of complex eigenvalues: ", len(complex_lambda))

#saving to csv
print("\nsaving to csv...")
df_real = pd.DataFrame(real_lambda, dtype = complex)
df_complex = pd.DataFrame(complex_lambda, dtype = complex)

df_complex.to_csv("Data/complex_eigenvalues_TL_N2.txt")
df_real.to_csv("Data/real_eigenvalues_TL_N2.txt")
