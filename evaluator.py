import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath


#Defining matricies
Lambda = pd.DataFrame([])
dif = pd.DataFrame([])
real_lambda_N2 = pd.DataFrame([])
complex_lambda_N2 = pd.DataFrame([])
zero_lambda = pd.DataFrame([])

#reading data
print("\nreading data")
#N=2
complex_lambda_N2 = pd.read_csv("Data/complex_eigenvalues_N2.txt", sep = ",", header = 1, \
    names = ["index", "number"])
complex_lambda_N2["number"] = complex_lambda_N2["number"].apply(lambda x: np.complex(x))
complex_lambda_N2 = complex_lambda_N2["number"].to_numpy()

real_lambda_N2 = pd.read_csv("Data/real_eigenvalues_N2.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N2["number"] = real_lambda_N2["number"].apply(lambda x: np.complex(x))
real_lambda_N2 = real_lambda_N2["number"].to_numpy()


#N=4
complex_lambda_N4 = pd.read_csv("Data/complex_eigenvalues_N4.txt", sep = ",", header = 1, \
    names = ["index", "number"])
complex_lambda_N4["number"] = complex_lambda_N4["number"].apply(lambda x: np.complex(x))
complex_lambda_N4 = complex_lambda_N4["number"].to_numpy()

real_lambda_N4 = pd.read_csv("Data/real_eigenvalues_N4.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N4["number"] = real_lambda_N4["number"].apply(lambda x: np.complex(x))
real_lambda_N4 = real_lambda_N4["number"].to_numpy()

#N=8
complex_lambda_N8 = pd.read_csv("Data/complex_eigenvalues_N8.txt", sep = ",", header = 1, \
    names = ["index", "number"])
complex_lambda_N8["number"] = complex_lambda_N8["number"].apply(lambda x: np.complex(x))
complex_lambda_N8 = complex_lambda_N8["number"].to_numpy()

real_lambda_N8 = pd.read_csv("Data/real_eigenvalues_N8.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N8["number"] = real_lambda_N8["number"].apply(lambda x: np.complex(x))
real_lambda_N8 = real_lambda_N8["number"].to_numpy()

#plotting
print("\nplotting...")
#N=2
#plotting histogram complex eigenvalues
plt.hist2d(np.append(complex_lambda_N2.real, real_lambda_N2.real), np.append(complex_lambda_N2.imag, \
real_lambda_N2.imag), bins = (201, 201), range = [[-4e5, -1e-6], [-1e-9, 1e-9]], density = True, cmap=plt.cm.BuPu)

plt.colorbar()
plt.title('Distr. compl. eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$Re \lambda$', fontsize = 13)
plt.xticks([-4e5, -3e5, -2e5, -1e5, 0])
plt.ylabel('$Im \lambda$', fontsize = 13)
plt.text(-3.8e5, 0.75e-9, '#compl. eigenvalues:   ' + str(len(np.append(complex_lambda_N2.real, real_lambda_N2.real))))
plt.text(-3.8e5, 0.6e-9, '$N^2$ = 16')
plt.savefig('Plots/Hist_N2.png', dpi=300)
plt.clf()

#histogram real eigenvalues
plt.hist(np.append(complex_lambda_N2.real, real_lambda_N2.real), range =[-4e5, -1e-6], bins = 1001, density = True)

plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda_{R}$', fontsize = 13)
plt.xticks([-4e5, -3e5, -2e5, -1e5, 0])
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-3.8e5, 4e-6, '#real eigenvalues:   ' + str(len(np.append(complex_lambda_N2.real, real_lambda_N2.real))))
plt.text(-3.8e5, 3.5e-6, '$N^2$ = 16')
plt.savefig('Plots/Hist_real_N2.png', dpi=300)
plt.clf()



#N=4
#plotting histogram complex eigenvalues
plt.hist2d(np.append(complex_lambda_N4.real, real_lambda_N4.real), np.append(complex_lambda_N4.imag, \
real_lambda_N4.imag), bins = (201, 201), range = [[-4e5, -1e-6], [-1e-9, 1e-9]], density = True, cmap=plt.cm.BuPu)

plt.colorbar()
plt.title('Distr. compl. eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$Re \lambda$', fontsize = 13)
plt.xticks([-4e5, -3e5, -2e5, -1e5, 0])
plt.ylabel('$Im \lambda$', fontsize = 13)
plt.text(-3.8e5, 0.75e-9, '#compl. eigenvalues:   ' + str(len(np.append(complex_lambda_N4.real, real_lambda_N4.real))))
plt.text(-3.8e5, 0.6e-9, '$N^2$ = 16')
plt.savefig('Plots/Hist_N4.png', dpi=300)
plt.clf()

#histogram real eigenvalues
plt.hist(np.append(complex_lambda_N4.real, real_lambda_N4.real), range =[-4e5, -1e-6], bins = 1001, density = True)

plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda_{R}$', fontsize = 13)
plt.xticks([-4e5, -3e5, -2e5, -1e5, 0])
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-3.8e5, 4e-6, '#real eigenvalues:   ' + str(len(np.append(complex_lambda_N4.real, real_lambda_N4.real))))
plt.text(-3.8e5, 3.5e-6, '$N^2$ = 16')
plt.savefig('Plots/Hist_real_N4.png', dpi=300)
plt.clf()

#N=8
#plotting histogram complex eigenvalues
plt.hist2d(np.append(complex_lambda_N8.real, real_lambda_N8.real), np.append(complex_lambda_N8.imag, \
real_lambda_N8.imag), bins = (201, 201), range = [[-4e5, -1e-6], [-1e-9, 1e-9]], density = True, cmap=plt.cm.BuPu)

plt.colorbar()
plt.title('Distr. compl. eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$Re \lambda$', fontsize = 13)
plt.xticks([-4e5, -3e5, -2e5, -1e5, 0])
plt.ylabel('$Im \lambda$', fontsize = 13)
plt.text(-3.8e5, 0.75e-9, '#compl. eigenvalues:   ' + str(len(np.append(complex_lambda_N8.real, real_lambda_N8.real))))
plt.text(-3.8e5, 0.6e-9, '$N^2$ = 16')
plt.savefig('Plots/Hist_N8.png', dpi=300)
plt.clf()

#histogram real eigenvalues
plt.hist(np.append(complex_lambda_N8.real, real_lambda_N8.real), range =[-4e5, -1e-6], bins = 1001, density = True)

plt.title('Distribution function real eigenvalues $P_{\lambda}$', fontsize = 15)
plt.xlabel('$\lambda_{R}$', fontsize = 13)
plt.xticks([-4e5, -3e5, -2e5, -1e5, 0])
plt.ylabel('probability distribution', fontsize = 13)
plt.text(-3.8e5, 4e-6, '#real eigenvalues:   ' + str(len(np.append(complex_lambda_N8.real, real_lambda_N8.real))))
plt.text(-3.8e5, 3.5e-6, '$N^2$ = 16')
plt.savefig('Plots/Hist_real_N8.png', dpi=300)
plt.clf()