import numpy as np
from scipy import linalg as lg


#Defining important constants

N = int(4)



#Defining matricies
S = np.zeros(shape=(N,N))
mathcalS = np.zeros(shape=(N^2,N^2))
mathcalA = np.zeros(shape=(N*N,N*N))


#Defining S and \mathcalS
for i in range(int(N/2)):
    S[i][i] = 1
for i in range(int(N/2), N):
    S[i][i] = -1

mathcalS = np.kron(S, S)


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

print(mathcalA)
