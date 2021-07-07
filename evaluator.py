import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath
from scipy import optimize


#---------------------------------reading data---------------------------------------
print("\nreading data...")
real_lambda_N3 = pd.read_csv("Data/real_eigenvalues_N3.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N3["number"] = real_lambda_N3["number"].apply(lambda x: np.complex(x))
real_lambda_N3 = real_lambda_N3["number"].to_numpy()
real_lambda_N4 = pd.read_csv("Data/real_eigenvalues_N4.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N4["number"] = real_lambda_N4["number"].apply(lambda x: np.complex(x))
real_lambda_N4 = real_lambda_N4["number"].to_numpy()
real_lambda_N5 = pd.read_csv("Data/real_eigenvalues_N5.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N5["number"] = real_lambda_N5["number"].apply(lambda x: np.complex(x))
real_lambda_N5 = real_lambda_N5["number"].to_numpy()
real_lambda_N6 = pd.read_csv("Data/real_eigenvalues_N6.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N6["number"] = real_lambda_N6["number"].apply(lambda x: np.complex(x))
real_lambda_N6 = real_lambda_N6["number"].to_numpy()
real_lambda_N7 = pd.read_csv("Data/real_eigenvalues_N7.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N7["number"] = real_lambda_N7["number"].apply(lambda x: np.complex(x))
real_lambda_N7 = real_lambda_N7["number"].to_numpy()
real_lambda_N8 = pd.read_csv("Data/real_eigenvalues_N8.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N8["number"] = real_lambda_N8["number"].apply(lambda x: np.complex(x))
real_lambda_N8 = real_lambda_N8["number"].to_numpy()
real_lambda_N9 = pd.read_csv("Data/real_eigenvalues_N9.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N9["number"] = real_lambda_N9["number"].apply(lambda x: np.complex(x))
real_lambda_N9 = real_lambda_N9["number"].to_numpy()


dif3 = pd.read_csv("Data/dif3.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif3 = dif3["number"].to_numpy()
dif4 = pd.read_csv("Data/dif4.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif4 = dif4["number"].to_numpy()
dif5 = pd.read_csv("Data/dif5.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif5 = dif5["number"].to_numpy()
dif6 = pd.read_csv("Data/dif6.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif6 = dif6["number"].to_numpy()
dif7 = pd.read_csv("Data/dif7.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif7 = dif7["number"].to_numpy()
dif8 = pd.read_csv("Data/dif8.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif8 = dif8["number"].to_numpy()
dif9 = pd.read_csv("Data/dif9.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif9 = dif9["number"].to_numpy()



#---------------------------------plotting eigenvalues-----------------------------
print("\nplotting...")
#------------------------------------N=3----------------------------------
plt.hist(real_lambda_N3.real, bins = 301, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-25, 0.175, '#real eigenvalues:   ' + str(len(real_lambda_N3.real)))
plt.text(-25, 0.165,  '$N^2$ = 9')
plt.savefig('Plots/Hist_real_N3.png', dpi=300)
plt.clf()
#-------------------------N=4--------------------------------------------
plt.hist(real_lambda_N4.real, bins = 301, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-30, 0.12, '#real eigenvalues:   ' + str(len(real_lambda_N4.real)))
plt.text(-30, 0.11, '$N^2$ = 16')
plt.savefig('Plots/Hist_real_N4.png', dpi=300)
plt.clf()

#------------------------------------------N=5---------------------------------------------
plt.hist(real_lambda_N5.real, bins = 301, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-50, 0.07, '#real eigenvalues:   ' + str(len(real_lambda_N5.real)))
plt.text(-50, 0.065, '$N^2$ = 25')
plt.savefig('Plots/Hist_real_N5.png', dpi=300)
plt.clf()

#--------------------------------------------------N=6----------------------------------
plt.hist(real_lambda_N6.real, bins = 301, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-60, 0.05, '#real eigenvalues:   ' + str(len(real_lambda_N6.real)))
plt.text(-60, 0.045,  '$N^2$ = 36')
plt.savefig('Plots/Hist_real_N6.png', dpi=300)
plt.clf()

#--------------------------------------------------N=7----------------------------------
plt.hist(real_lambda_N7.real, bins = 301, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-90, 0.03, '#real eigenvalues:   ' + str(len(real_lambda_N7.real)))
plt.text(-90, 0.028,  '$N^2$ = 49')
plt.savefig('Plots/Hist_real_N7.png', dpi=300)
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
plt.hist(real_lambda_N8.real, bins = 301, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-110, 0.03, '#real eigenvalues:   ' + str(len(real_lambda_N8.real)))
plt.text(-110, 0.028,  '$N^2$ = 64')
plt.savefig('Plots/Hist_real_N8.png', dpi=300)
plt.clf()

#--------------------------------------------------N=9----------------------------------
plt.hist(real_lambda_N9.real, bins = 301, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-150, 0.0225, '#real eigenvalues:   ' + str(len(real_lambda_N9.real)))
plt.text(-150, 0.02,  '$N^2$ = 81')
plt.savefig('Plots/Hist_real_N9.png', dpi=300)
plt.clf()



#-----------------------------plotting correlation-----------------------------
print("\nplotting correlation...")
#N=3
plt.hist(dif3, bins = 301, range = [0, 10], density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(6, 0.55, '#eigenvalues:   ' + str(len(dif3)))
plt.text(6, 0.5,  '$N^2$ = 9')
plt.savefig('Plots/correlation_N3.png', dpi=300)
plt.clf()

#N=4
plt.hist(dif4, bins = 301, range = [1e-12, 10], density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$ without zeros', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(6, 0.55, '#eigenvalues:   ' + str(len(dif4)))
plt.text(6, 0.5,  '$N^2$ = 16')
plt.savefig('Plots/correlation_N4.png', dpi=300)
plt.clf()

#N=5
plt.hist(dif5, bins = 301, density = True)
plt.errorbar(np.linspace(0,0.3,1000), 5.5*np.linspace(0,0.3,1000), color = "black", fmt = "--")
plt.title('Correlation eigenvalues $\Delta\lambda$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(10, 0.8, '#eigenvalues:   ' + str(len(dif5)))
plt.text(10, 0.7,  '$N^2$ = 25')
plt.legend(['generated eigenvalues', "linear fit: y = 5.5x"], fontsize = 13)
plt.savefig('Plots/correlation_N5.png', dpi=300)
plt.clf()

#N=6
plt.hist(dif6, bins = 301, density = True, range = [0,3])
plt.errorbar(np.linspace(0,0.3,1000), 5.5*np.linspace(0,0.3,1000), color = "black", fmt = "--")
plt.title('Correlation eigenvalues $\Delta\lambda$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(2, 0.8, '#eigenvalues:   ' + str(len(dif6)))
plt.text(2, 0.7,  '$N^2$ = 36')
plt.legend(['generated eigenvalues', "linear fit: y = 5.5x"], fontsize = 13)
plt.savefig('Plots/correlation_N6.png', dpi=300)
plt.clf()

#N=7
plt.hist(dif7, bins = 301, density = True, range = [0,3])
plt.errorbar(np.linspace(0,0.3,1000), 5.5*np.linspace(0,0.3,1000), color = "black", fmt = "--")
plt.title('Correlation eigenvalues $\Delta\lambda$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(0.7, 1.2, '#eigenvalues:   ' + str(len(dif7)))
plt.text(0.7, 1.4,  '$N^2$ = 49')
plt.legend(['generated eigenvalues', "linear fit: y = 5.5x"], fontsize = 13)
plt.savefig('Plots/correlation_N7.png', dpi=300)
plt.clf()

#N=8
plt.hist(dif8, bins = 301, range = [1e-12, 10], density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$ without zeros', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(5, 0.6, '#eigenvalues:   ' + str(len(dif8)))
plt.text(5, 0.55,  '$N^2$ = 64')
plt.savefig('Plots/correlation_N8.png', dpi=300)
plt.clf()

#N=9
plt.hist(dif9, bins = 301, range = [0, 10], density = True)
plt.errorbar(np.linspace(0,0.25,1000), 5.5*np.linspace(0,0.25,1000), color = "black", fmt = "--")
plt.title('Correlation eigenvalues $\Delta\lambda$ without zeros', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.text(5, 0.6, '#eigenvalues:   ' + str(len(dif9)))
plt.text(5, 0.55,  '$N^2$ = 81')
plt.legend(['generated eigenvalues', "linear fit: y = 5.5x"], fontsize = 13)
plt.savefig('Plots/correlation_N9.png', dpi=300)
plt.clf()


#-----------------------------------------mean lambda(N)-------------------------------------------
def exp(x, a, b):
    return a*np.e**(b*x)
def lin(x, b, c):
    return b*x**2 + c*x

mean_lambda = pd.DataFrame([float(7.0), float(17.0), float(27.0), float(45.0), float(65.0), 95], columns = ["value"])
N = pd.DataFrame([4, 5, 6, 7, 8, 9], columns = ["value"])
params, params_cov = optimize.curve_fit(lin, N["value"], mean_lambda["value"])
print(params)
plt.errorbar(N, mean_lambda, fmt = "D")
plt.errorbar(np.linspace(4,9, 1000), lin(np.linspace(4,9, 1000), params[0], params[1]), fmt = "--", color = "black")
plt.title(r'$\langle\lambda\rangle$ as function of N', fontsize = 15)
plt.xlabel('$N$', fontsize = 13)
plt.ylabel(r"$\langle \lambda \rangle$", fontsize = 13)
plt.legend(['generated eigenvalues', "polynomial fit: $y = 1.8x^2-6.0x$"], fontsize = 13)
plt.savefig('Plots/mean_lambda.png', dpi=300)
plt.clf()
