import numpy as np
from scipy import linalg as lg
import matplotlib.pyplot as plt
import pandas as pd
import cmath
from scipy import optimize


plt.style.use('bmh')

plt.style.use('Solarize_Light2')

plt.style.use('classic')

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
real_lambda_N10 = pd.read_csv("Data/real_eigenvalues_N10.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N10["number"] = real_lambda_N10["number"].apply(lambda x: np.complex(x))
real_lambda_N10 = real_lambda_N10["number"].to_numpy()
real_lambda_N16 = pd.read_csv("Data/real_eigenvalues_N16.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N16["number"] = real_lambda_N16["number"].apply(lambda x: np.complex(x))
real_lambda_N16 = real_lambda_N16["number"].to_numpy()
real_lambda_N20 = pd.read_csv("Data/real_eigenvalues_N20.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N20["number"] = real_lambda_N20["number"].apply(lambda x: np.complex(x))
real_lambda_N20 = real_lambda_N20["number"].to_numpy()
real_lambda_N48 = pd.read_csv("Data/real_eigenvalues_N48.txt", sep = ",", header = 0, \
    names = ["index", "number"])
real_lambda_N48["number"] = real_lambda_N48["number"].apply(lambda x: np.complex(x))
real_lambda_N48 = real_lambda_N48["number"].to_numpy()


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
dif10 = pd.read_csv("Data/dif10.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif10 = dif10["number"].to_numpy()
dif20 = pd.read_csv("Data/dif20.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif20 = dif20["number"].to_numpy()
dif48 = pd.read_csv("Data/dif48.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif48 = dif48["number"].to_numpy()





#all the same
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}

plt.rc('font', **font)
plt.gcf().subplots_adjust(bottom=0.15)
n,x = np.histogram(dif3, bins = 100, density = True, range = [0, 5])
plt.plot(x[:-1], n)
n,x= np.histogram(dif4, bins = 100, density = True, range = [0, 5])
plt.plot(x[:-1], n)
n,x= np.histogram(dif6, bins = 100, density = True, range = [0, 5])
plt.plot(x[:-1], n)
n,x= np.histogram(dif9, bins = 100, density = True, range = [0, 5])
plt.plot(x[:-1], n)
n,x= np.histogram(dif20, bins = 100, density = True, range = [0, 5])
plt.plot(x[:-1], n)
#plt.title('Eigenvalue correlation for different dimensions\n')
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$\n', fontsize=18)
plt.ylabel('probability distribution')
plt.legend(["N=3", "N=4", "N=6", "N=9", "N=20"])
plt.savefig('Plots/correlation_same.png', dpi=300)
plt.clf()

#Wigner
def wigner(s):
    return s/2 * np.exp(-s**2/4)
def poisson(s):
    return 0.5*np.exp(-s)

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}

plt.rc('font', **font)
plt.gcf().subplots_adjust(bottom=0.15)
plt.plot(np.linspace(0, 5, 1000), wigner(np.linspace(0, 5, 1000)))
plt.plot(np.linspace(0, 5, 1000), poisson(np.linspace(0, 5, 1000)), color = "black", linestyle = "--")
#plt.title('Eigenvalue correlation for different dimensions\n')
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize=18)
plt.ylabel('probability distribution')
plt.yticks(np.linspace(0,0.4,5))
plt.legend(["Wigner surmise", "poisson distribution"])
plt.savefig('Plots/wigner_surmise.png', dpi=300)
plt.clf()

#-----------------------------------------mean lambda(N)-------------------------------------------
def exp(x, a, b):
    return a*np.e**(b*x)
def lin(x, b, c):
    return b*x**c

mean_lambda = pd.DataFrame([6, 13, 28, 50, 82, 123, 175, 250, 1020, 2000], columns = ["value"])
N = pd.DataFrame([3, 4, 5, 6, 7, 8, 9, 10, 16, 20], columns = ["value"])
params, params_cov = optimize.curve_fit(lin, N["value"], mean_lambda["value"])
print(params)
plt.errorbar(N, mean_lambda, fmt = "D")
plt.errorbar(np.linspace(3,20, 1000), lin(np.linspace(3,20, 1000), params[0], params[1]), fmt = "--", color = "black")
plt.title(r'$\langle\lambda\rangle$ as function of N', fontsize = 15)
plt.xlabel('$N$', fontsize = 13)
plt.ylabel(r"$\langle \lambda \rangle$", fontsize = 13)
plt.legend(['generated eigenvalues', "polynomial fit: $y = 0.2x^2-3.0x$"], fontsize = 13)
plt.savefig('Plots/mean_lambda.png', dpi=300)
plt.clf()

#-------------------------stacked subplots-----------------------------
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

plt.rc('font', **font)
plt.rcParams["figure.figsize"] = (16,18)
fig, axs = plt.subplots(3,2)
axs[0,0].hist(real_lambda_N3.real, bins = 201, density = True)
axs[0,0].set_title('N=3')
axs[0,0].set_xticks([-30,-20,-10,0])
axs[0,0].set_yticks(np.linspace(0,0.12,4))
axs[0,0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[0,0].set_ylabel("probability distribution")
axs[0,1].hist(real_lambda_N4.real, bins = 201, density = True)
axs[0,1].set_title('N=4')
axs[0,1].set_xticks([-60,-40,-20,0])
axs[0,1].set_yticks(np.linspace(0,0.06,4))
axs[0,1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[1,0].hist(real_lambda_N6.real, bins = 201, density = True)
axs[1,0].set_title('N=6')
axs[1,0].set_xticks(np.linspace(-100,-20,5))
axs[1,0].set_yticks(np.linspace(0,0.03,4))
axs[1,0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[1,0].set_ylabel("probability distribution\n")
axs[1,1].hist(real_lambda_N9.real, bins = 201, density = True)
axs[1,1].set_title('N=9')
axs[1,1].set_xticks(np.linspace(-270,-120,4))
axs[1,1].set_yticks(np.linspace(0,0.012,4))
axs[1,1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[2,0].hist(real_lambda_N10.real, bins = 201, density = True, range = [-360,-180])
axs[2,0].set_title('N=10')
axs[2,0].set_xticks(np.linspace(-330,-180,4))
axs[2,0].set_yticks(np.linspace(0,0.012,4))
axs[2,0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[2,0].set_xlabel("$\lambda$")
axs[2,0].set_ylabel("probability distribution")
axs[2,1].hist(real_lambda_N20.real, bins = 201, density = True)
axs[2,1].set_title('N=20')
axs[2,1].set_xticks(np.linspace(-2300,-1700,4))
axs[2,1].set_yticks(np.linspace(0,0.003,4))
axs[2,1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[2,1].set_xlabel("$\lambda$")
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig('Plots/test.png', dpi=300)
plt.clf()
