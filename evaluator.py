import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath



#---------------------------------reading data---------------------------------------
print("\nreading data...")
#N=4
real_lambda_N4 = pd.read_csv("Data/real_eigenvalues_N4.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N4["number"] = real_lambda_N4["number"].apply(lambda x: np.complex(x))
real_lambda_N4 = real_lambda_N4["number"].to_numpy()
real_lambda_N5 = pd.read_csv("Data/real_eigenvalues_N5.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N5["number"] = real_lambda_N5["number"].apply(lambda x: np.complex(x))
real_lambda_N5 = real_lambda_N5["number"].to_numpy()
real_lambda_N8 = pd.read_csv("Data/real_eigenvalues_N8.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N8["number"] = real_lambda_N8["number"].apply(lambda x: np.complex(x))
real_lambda_N8 = real_lambda_N8["number"].to_numpy()


Lambda_sort_N4 = pd.read_csv("Data/sorted_N4.txt", sep = ",", header = 1, \
    names = ["index", "number"]).to_numpy()
Lambda_sort_N5 = pd.read_csv("Data/sorted_N5.txt", sep = ",", header = 1, \
    names = ["index", "number"]).to_numpy()
Lambda_sort_N8 = pd.read_csv("Data/sorted_N8.txt", sep = ",", header = 1, \
    names = ["index", "number"]).to_numpy()

dif4 = pd.read_csv("Data/dif4.txt", sep = ",", header = 1, \
    names = ["index", "number"])
dif4 = dif4["number"].to_numpy()
#dif5 = pd.read_csv("Data/dif5.txt").to_numpy()
#dif8 = pd.read_csv("Data/dif8.txt").to_numpy()



#---------------------------------plotting-----------------------------
print("\nplotting...")
#-------------------------N=4--------------------------------------------
#plotting histogram complex eigenvalues
plt.hist2d(real_lambda_N4.real, real_lambda_N4.imag, bins = (201, 201), range = [[-20, 0], \
[-1e-12, 1e-12]], density = True, cmap=plt.cm.BuPu)
plt.colorbar()
plt.title('Distr. compl. eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$Re \lambda$', fontsize = 13)
plt.ylabel('$Im \lambda$', fontsize = 13)
plt.text(-18, 0.75e-12, '#compl. eigenvalues:   ' + str(len(real_lambda_N4.real)))
plt.text(-18, 0.6e-12, '$N^2$ = 16')
plt.savefig('Plots/Hist_N4.png', dpi=300)
plt.clf()

#histogram real eigenvalues
plt.hist(real_lambda_N4.real, bins = 501, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda_{R}$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-30, 0.12, '#real eigenvalues:   ' + str(len(real_lambda_N4.real)))
plt.text(-30, 0.11, '$N^2$ = 16')
plt.savefig('Plots/Hist_real_N4.png', dpi=300)
plt.clf()

#------------------------------------------N=5---------------------------------------------
#plotting histogram complex eigenvalues
plt.hist2d(real_lambda_N5.real, real_lambda_N5.imag, bins = (201, 201), \
range = [[-45, 0], [-1e-12, 1e-12]], density = True, cmap=plt.cm.BuPu)
plt.colorbar()
plt.title('Distr. compl. eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$Re \lambda$', fontsize = 13)
plt.ylabel('$Im \lambda$', fontsize = 13)
plt.text(-42, 0.75e-12, '#compl. eigenvalues:   ' + str(len(real_lambda_N5.real)))
plt.text(-42, 0.6e-12, '$N^2$ = 16')
plt.savefig('Plots/Hist_N5.png', dpi=300)
plt.clf()

#histogram real eigenvalues
plt.hist(real_lambda_N5.real, bins = 501, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda_{R}$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-50, 0.07, '#real eigenvalues:   ' + str(len(real_lambda_N5.real)))
plt.text(-50, 0.065, '$N^2$ = 16')
plt.savefig('Plots/Hist_real_N5.png', dpi=300)
plt.clf()


#----------------------------------------------N=8-----------------------------------------
#plotting histogram complex eigenvalues
plt.hist2d(real_lambda_N8.real,real_lambda_N8.imag, bins = (201, 201), density = True, \
range = [[-130, 0], [-1e-12, 1e-12]], cmap=plt.cm.BuPu)
plt.colorbar()
plt.title('Distr. compl. eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$Re \lambda$', fontsize = 13)
plt.ylabel('$Im \lambda$', fontsize = 13)
plt.text(-120, 0.75e-12, '#compl. eigenvalues:   ' + str(len(real_lambda_N8.real)))
plt.text(-120, 0.65e-12, '$N^2$ = 64')
plt.savefig('Plots/Hist_N8.png', dpi=300)
plt.clf()

#histogram real eigenvalues
plt.hist(real_lambda_N8.real, bins = 501, density = True)

plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda_{R}$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-120, 0.03, '#real eigenvalues:   ' + str(len(real_lambda_N8.real)))
plt.text(-120, 0.028,  '$N^2$ = 64')
plt.savefig('Plots/Hist_real_N8.png', dpi=300)
plt.clf()



#-----------------------------plotting correlation-----------------------------
print("\nplotting correlation...")
#N=4
plt.hist(dif4, bins = 301, range = [1e-5, 10], density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-120, 0.03, '#real eigenvalues:   ' + str(len(dif4)))
plt.text(-120, 0.028,  '$N^2$ = 16')
plt.savefig('Plots/correlation_N4.png', dpi=300)
plt.clf()

#N=5
plt.hist(dif5, bins = 301, density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-120, 0.03, '#real eigenvalues:   ' + str(len(dif5)))
plt.text(-120, 0.028,  '$N^2$ = 25')
plt.savefig('Plots/correlation_N5.png', dpi=300)
plt.clf()
'''
#N=8
plt.hist(dif8, bins = 501, density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-120, 0.03, '#real eigenvalues:   ' + str(len(dif8)))
plt.text(-120, 0.028,  '$N^2$ = 64')
plt.savefig('Plots/correlation_N8.png', dpi=300)
plt.clf()
'''
