import numpy as np
from scipy import linalg as lg
from scipy import signal
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
dif16 = pd.read_csv("Data/dif16.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif16 = dif16["number"].to_numpy()
dif20 = pd.read_csv("Data/dif20.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif20 = dif20["number"].to_numpy()





#---------------------------------plotting eigenvalues-----------------------------
print("\nplotting...")
#------------------------------------N=3----------------------------------
plt.hist(real_lambda_N3.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=3$', fontsize = 15)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N3.real))], loc = 2)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.savefig('Plots/Hist_real_N3.png', dpi=300)
plt.clf()
#-------------------------N=4--------------------------------------------
plt.hist(real_lambda_N4.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=4$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N4.real))], loc = 2)
plt.savefig('Plots/Hist_real_N4.png', dpi=300)
plt.clf()

#------------------------------------------N=5---------------------------------------------
plt.hist(real_lambda_N5.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=5$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N5.real))], loc = 2)
plt.ylabel('probability distribution', fontsize = 13)
plt.savefig('Plots/Hist_real_N5.png', dpi=300)
plt.clf()

#--------------------------------------------------N=6----------------------------------
plt.hist(real_lambda_N6.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=6$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N6.real))], loc = 2)
plt.savefig('Plots/Hist_real_N6.png', dpi=300)
plt.clf()

#--------------------------------------------------N=7----------------------------------
plt.hist(real_lambda_N7.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, , $N=7$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N7.real))], loc = 2)
plt.savefig('Plots/Hist_real_N7.png', dpi=300)
plt.clf()


#----------------------------------------------N=8-----------------------------------------
plt.hist(real_lambda_N8.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=8$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N8.real))], loc = 2)
plt.savefig('Plots/Hist_real_N8.png', dpi=300)
plt.clf()

#--------------------------------------------------N=9----------------------------------
plt.hist(real_lambda_N9.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=9$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N9.real))], loc = 2)
plt.savefig('Plots/Hist_real_N9.png', dpi=300)
plt.clf()

#--------------------------------------------------N=10----------------------------------
plt.hist(real_lambda_N10.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=10$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N10.real))], loc = 2)
plt.savefig('Plots/Hist_real_N10.png', dpi=300)
plt.clf()

#--------------------------------------------------N=16----------------------------------
plt.hist(real_lambda_N16.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=16$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N16.real))], loc = 2)
plt.savefig('Plots/Hist_real_N16.png', dpi=300)
plt.clf()

#--------------------------------------------------N=20----------------------------------
plt.hist(real_lambda_N20.real, bins = 201, density = True)
plt.title('Distribution function real eigenvalues $P_{\lambda}$, $N=20$', fontsize = 15)
plt.xlabel('$\lambda$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated eigenvalues: "+ str(len(real_lambda_N20.real))], loc = 2)
plt.savefig('Plots/Hist_real_N20.png', dpi=300)
plt.clf()



#-----------------------------plotting correlation-----------------------------
print("\nplotting correlation...")
#N=3
plt.hist(dif3, bins = 201, density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=3$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif3))], loc = 1)
plt.savefig('Plots/correlation_N3.png', dpi=300)
plt.clf()

#N=4
plt.hist(dif4, bins = 201, density = True)#, range = [0, 3])
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=4$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif4))], loc = 1)
plt.savefig('Plots/correlation_N4.png', dpi=300)
plt.clf()

#N=5
plt.hist(dif5, bins = 201, density = True, range = [0,5])
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=5$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif5))], loc = 1)
plt.savefig('Plots/correlation_N5.png', dpi=300)
plt.clf()

#N=6
plt.hist(dif6, bins = 201, density = True, range = [0,5])
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=6$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif6))], loc = 1)
plt.savefig('Plots/correlation_N6.png', dpi=300)
plt.clf()

#N=7
plt.hist(dif7, bins = 201, density = True, range = [0,5])
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=7$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif7))], loc = 1)
plt.savefig('Plots/correlation_N7.png', dpi=300)
plt.clf()

#N=8
plt.hist(dif8, bins = 201, range = [0, 10], density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=8$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif8))], loc = 1)
plt.savefig('Plots/correlation_N8.png', dpi=300)
plt.clf()

#N=9
plt.hist(dif9, bins = 201, range = [0, 10], density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=9$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif9))], loc = 1)
plt.savefig('Plots/correlation_N9.png', dpi=300)
plt.clf()

#N=10
plt.hist(dif10, bins = 201, density = True, range = [0,5])
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=10$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif10))], loc = 1)
plt.savefig('Plots/correlation_N10.png', dpi=300)
plt.clf()

#N=16
plt.hist(dif16, bins = 201, range = [0, 10], density = True)
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=16$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif16))], loc = 1)
plt.savefig('Plots/correlation_N16.png', dpi=300)
plt.clf()

#N=20
plt.hist(dif20, bins = 201, density = True, range = [0,5])
plt.title('Correlation eigenvalues $\Delta\lambda$, $N=20$', fontsize = 15)
plt.xlabel('$|\lambda_i - \lambda_{i+1}|$', fontsize = 13)
plt.ylabel('probability distribution', fontsize = 13)
plt.legend(["Calculated nearest neighbour distances: "+ str(len(dif20))], loc = 1)
plt.savefig('Plots/correlation_N20.png', dpi=300)
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
