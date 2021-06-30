import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath


#Defining matricies
Lambda = pd.DataFrame([])
dif = pd.DataFrame([])
real_lambda = pd.DataFrame([])
complex_lambda_N2 = pd.DataFrame([])
zero_lambda = pd.DataFrame([])

#reading data
#real_lambda = np.loadtxt("Data/real_eigenvalues.txt").view(complex)
complex_lambda_N2 = pd.read_csv("Data/complex_eigenvalues_N2.txt", sep = ",", header = 1, \
    names = ["index", "number"])
complex_lambda_N2["number"] = complex_lambda_N2["number"].apply(lambda x: np.complex(x))
complex_lambda_N2 = complex_lambda_N2["number"].to_numpy()

real_lambda = pd.read_csv("Data/real_eigenvalues_N2.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda["number"] = real_lambda["number"].apply(lambda x: np.complex(x))
real_lambda = real_lambda["number"].to_numpy()


#plotting histogram
plt.hist2d(complex_lambda_N2.real, complex_lambda_N2.imag, bins = (301, 301), range = [[-1000, 1000],[-1000, 1000]], density = True, cmap=plt.cm.BuPu)

plt.colorbar()
plt.title('Distribution fuction $P_{\lambda}$', fontsize = 15)
plt.xlabel('$Re \lambda$', fontsize = 13)
plt.ylabel('$Im \lambda$', fontsize = 13)
plt.text(-850, 850, '#compl. eigenvalues:   ' + str(len(complex_lambda_N2.real)))
plt.text(-850, 750, '$N^2$ = 4')
plt.savefig('Plots/Hist.png', dpi=300)
plt.clf()

'''
#distance
plt.hist(dist, bins = 601, range = [0, 0.04], density = True)

plt.title('Histogram correlation', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.savefig('Plots/Hist_dist.png', dpi=300)
plt.clf()
'''
