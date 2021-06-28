import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd


#Defining important constants

N = int(10)
M = int(80)


#Defining matricies
mathcalA = np.zeros(shape=(N*N,N*N))
kappa = np.zeros(shape=(N*N,N*N))
U = np.zeros(shape=(N*N,N*N))
tildekappa = np.zeros(shape=(N*N,N*N))
K_0 = np.zeros(shape=(N*N,N*N))
kappa_0 = np.zeros(shape=(N*N,N*N))
mathcalL = np.zeros(shape=(N*N,N*N))
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

    #Calculating mathcalL from K_0 via Lange,Timm (easy way, by back-transforming) Eq.:(23)
    kappa_0 = U.T.dot(K_0.dot(U))
    for i in range(N):
        for j in range(N):
            for m in range(N):
                for n in range(N):
                    for sum in range(N):
                        mathcalL[N*i+j][N*m+n] += signal.unit_impulse(N*N, i)*0.5*kappa_0[sum*N+j][sum*N+n] + signal.unit_impulse(N*N, j)*0.5*kappa_0[i*N+sum][m*N+sum]
                    mathcalL[N*i+j][N*m+n] += kappa_0[N*i+m][N*j+n]



    values = lg.eigvals(mathcalL)

    Lambda = np.append(Lambda, values)


print(Lambda, len(Lambda))
Lambda_sort = np.sort(Lambda)

for i in range(len(Lambda_sort)-1):
     dif = np.append(dif, Lambda_sort[i+1]-Lambda_sort[i])

print("Das ist die Differenz:   ", dif)

#plotting histogram
n, bins, patches = plt.hist(Lambda, bins = 1000, range = (-1, 10), density = True)

plt.title('Histogram eigenvalues', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.savefig('Plots/Hist.png', dpi=300)
plt.clf()

#distance
n, bins, patches = plt.hist(dif, bins = 1000, range = (0, 0.0002), density = True)

plt.title('Histogram correlation', fontsize = 15)
plt.xlabel('$\lambda_i - \lambda_j$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.savefig('Plots/Hist_dist.png', dpi=300)
plt.clf()
