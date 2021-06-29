import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath


#Defining important constants

N = int(4)
M = int(8000)


#Defining matricies
mathcalA = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
kappa = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
U = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
tildekappa = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
K_0 = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
kappa_0 = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
mathcalL = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
Lambda = np.array([])
dif = np.array([])

A = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
B = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
C = np.zeros(shape=(N*N,N*N),dtype=np.complex_)
D = np.zeros(shape=(N*N,N*N),dtype=np.complex_)


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

    #kappa = AA^T
    kappa = mathcalA.dot(mathcalA.T)

    #tildekappa = UkappaU^T
    tildekappa = U.dot(kappa.dot(U.T))

    K_0 = tildekappa
    for a in range(N*N):
            K_0[N*(N-1)+(N-1)][a] = 0
    for b in range(N*N):
            K_0[b][N*(N-1)+(N-1)] = 0

    '''#Calculating mathcalL from K_0 via Lange,Timm (easy way, by back-transforming) Eq.:(23)
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
                    mathcalL[N*i+j][N*m+n] += kappa_0[N*i+m][N*j+n]'''

    #Calculating mathcalL using Timm, Lange Appendix (A1)-(A5)
    #A
    for m in range(N):
        for n in range(N):
            for p in range(N):
                for q in range(N):

                    if m == p:
                        for sum in range(N):
                            if sum != n:
                                if sum != q:
                                    A[N*m+n][N*p+q] += 0.5*K_0[sum*N+n][sum*N+q]
                    if n == p:
                        for sum in range(N):
                            if sum != m:
                                if sum != p:
                                    A[N*m+n][N*p+q] += 0.5*K_0[sum*N+p][sum*N+m]
                    if m != p:
                        if n!= q:
                            A[N*m+n][N*p+q] += K_0[m*N+p][n*N+m]
    #B
    for m in range(N):
        for n in range(N):
            for p in range(N):
                for q in range(N):

                    if m != p:
                        if n == q:
                            B[N*m+n][N*p+q] += (-np.sqrt((n)/(n+1))*K_0[N*m+p][N*(n-1)+(n-1)])
                            for k in range(n, N-1):
                                B[N*m+n][N*p+q] += K_0[m*N+p][k*N+k]/np.sqrt((k+1)*(k+2))

                    if m == p:
                        if n != q:
                            B[N*m+n][N*p+q] += 0.5*(np.sqrt((q)/(q+1))*K_0[N*q+n][N*(q-1)+(q-1)])
                            for k in range(q, N-1):
                                B[N*m+n][N*p+q] += (-0.5*K_0[q*N+n][k*N+k]/np.sqrt((k+1)*(k+2)))

                    if m != p:
                        if n == q:
                            B[N*m+n][N*p+q] += 0.5*(np.sqrt((m)/(m+1))*K_0[N*m+p][N*(m-1)+(m-1)])
                            for k in range(m, N-1):
                                B[N*m+n][N*p+q] += (-0.5*K_0[m*N+p][k*N+k]/np.sqrt((k+1)*(k+2)))

    #C
    for m in range(N):
        for n in range(N):
            for p in range(N):
                for q in range(N):

                    if m == p:
                        if n != q:
                            C[N*m+n][N*p+q] += (-np.sqrt((m)/(m+1))*K_0[N*m+p][N*(m-1)+(m-1)])
                            for k in range(m, N-1):
                                C[N*m+n][N*p+q] += K_0[k*N+k][n*N+q]/np.sqrt((k+1)*(k+2))

                    if m == p:
                        if n != q:
                            C[N*m+n][N*p+q] += 0.5*(np.sqrt((p)/(p+1))*K_0[N*(p-1)+(p-1)][N*p+m])
                            for k in range(p, N-1):
                                C[N*m+n][N*p+q] += (-0.5*K_0[k*N+k][p*N+m]/np.sqrt((k+1)*(k+2)))

                    if m == p:
                        if n != q:
                            C[N*m+n][N*p+q] += 0.5*(np.sqrt((n)/(n+1))*K_0[N*(n-1)+(n-1)][N*n+q])
                            for k in range(n, N-1):
                                C[N*m+n][N*p+q] += (-0.5*K_0[k*N+k][n*N+q]/np.sqrt((k+1)*(k+2)))
    #D
    for m in range(N):
        for n in range(N):
            for p in range(N):
                for q in range(N):

                    if m == p:
                        if n == q:
                            D[N*m+n][N*p+q] += (np.sqrt((m*n)/((m+1)*(n+1)))*K_0[N*(m-1)+(m-1)][N*(n-1)+(n-1)])
                            for l in range(n, N-1):
                                D[N*m+n][N*p+q] += -np.sqrt(m/(m+1))*K_0[(m-1)*N+(m-1)][l*N+l]/np.sqrt((l+1)*(l+2))

                            for k in range(m, N-1):
                                D[N*m+n][N*p+q] += -np.sqrt(n/(n+1))*K_0[k*N+k][(n-1)*N+(n-1)]/np.sqrt((k+1)*(k+2))

                            for k in range(m, N-1):
                                for l in range(n, N-1):
                                    D[N*m+n][N*p+q] += K_0[N*k+k][N*l+l]/np.sqrt((k+1)*(k+2)*(l+1)*(l+2))

                    D[N*m+n][N*p+q] += -m/(m+1) * K_0[N*(m-1)+(m-1)][N*(m-1)+(m-1)]
                    for k in range(m, N-1):
                        D[N*m+n][N*p+q] += np.sqrt(m/(m+1))*(K_0[(m-1)*N+(m-1)][k*N+k] + K_0[k*N+k][(m-1)*N+(m-1)])/np.sqrt((k+1)*(k+2))
                    for k in range(m, N-1):
                        for l in range(m, N-1):
                            D[N*m+n][N*p+q] += -K_0[N*k+k][N*l+l]/np.sqrt((k+1)*(k+2)*(l+1)*(l+2))

                    D[N*m+n][N*p+q] += -n/(n+1) * K_0[N*(n-1)+(n-1)][N*(n-1)+(n-1)]
                    for k in range(n, N-1):
                        D[N*m+n][N*p+q] += np.sqrt(n/(n+1))*(K_0[(n-1)*N+(n-1)][k*N+k] + K_0[k*N+k][(n-1)*N+(n-1)])/np.sqrt((k+1)*(k+2))
                    for k in range(n, N-1):
                        for l in range(n, N-1):
                            D[N*m+n][N*p+q] += -K_0[N*k+k][N*l+l]/np.sqrt((k+1)*(k+2)*(l+1)*(l+2))






    mathcalL = A+B+C+D
    values = lg.eigvals(mathcalL)

    Lambda = np.append(Lambda, values)


print(Lambda, len(Lambda))
Lambda_sort = np.sort(Lambda)

Lambda_real = Lambda.real
Lambda_imag = Lambda.imag

for i in range(len(Lambda_sort)-1):
     dif = np.append(dif, Lambda_sort[i+1]-Lambda_sort[i])
     dist= np.sqrt(dif.real**2 + dif.imag**2)


print("Das ist die Differenz:   ", dif)

#plotting histogram
plt.hist2d(Lambda_real, Lambda_imag, bins = (301, 301), range = [[-1000, 1000],[-1000, 1000]], density = True)

plt.colorbar()
plt.title('Histogram eigenvalues', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.savefig('Plots/Hist.png', dpi=300)
plt.clf()

#distance
plt.hist(dist, bins = 1000, density = True)

plt.title('Histogram correlation', fontsize = 15)
plt.xlabel('$\lambda_i - \lambda_j$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.savefig('Plots/Hist_dist.png', dpi=300)
plt.clf()
