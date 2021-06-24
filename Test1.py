import numpy as np
from scipy import linalg as lg


#Defining a matrix A, B

A = np.array([  [1,2] , [3,4]   ])

B = np.array([  [6,1] , [5,1]   ])

#Carrying out calculations

sum = A+B

dif = A-B

prod = A.dot(B)         #Matrix multiplication

transpose = A.T

determinantB = lg.det(B)

inverseB = lg.inv(B)

values, vectors = lg.eig(B)

print(A)

print(B)
