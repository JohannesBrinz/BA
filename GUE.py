import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath
from skrmt.ensemble import GaussianEnsemble

values1 = []
values1 = []
Lambda = []
dif = []

M = int(5e3)
N = int(20)

for i in range(M):
    GOE1 = GaussianEnsemble(beta=1, n=int(N**2/2))
    GOE2 = GaussianEnsemble(beta=1, n=int(N**2/2))

    values1 = lg.eigh(GOE1.matrix, eigvals_only=True)
    values2 = lg.eigh(GOE2.matrix, eigvals_only=True)
    Lambda = np.append(Lambda, values1)
    Lambda = np.append(Lambda, values2)
    values = np.append(values1, values2)
    values = np.sort(values)
    for i in range(len(values)-2):
         dif = np.append(dif, values[i+1]-values[i])

plt.hist(dif, bins=201, range=[0,2], density=True)
plt.savefig('Plots/GOE.png', dpi=300)
plt.clf()

df_dif = pd.DataFrame(dif.real)
df_dif.to_csv("Data/dif_GOE_20.txt")
