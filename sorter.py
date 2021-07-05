import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath



#reading data
print("\nreading data...")
#N=4
real_lambda_N4 = pd.read_csv("Data/real_eigenvalues_N4.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N4["number"] = real_lambda_N4["number"].apply(lambda x: np.complex(x))
real_lambda_N4 = real_lambda_N4["number"].to_numpy()

#N=5
real_lambda_N5 = pd.read_csv("Data/real_eigenvalues_N5.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N5["number"] = real_lambda_N5["number"].apply(lambda x: np.complex(x))
real_lambda_N5 = real_lambda_N5["number"].to_numpy()

#N=8
real_lambda_N8 = pd.read_csv("Data/real_eigenvalues_N8.txt", sep = ",", header = 1, \
    names = ["index", "number"])
real_lambda_N8["number"] = real_lambda_N8["number"].apply(lambda x: np.complex(x))
real_lambda_N8 = real_lambda_N8["number"].to_numpy()


#sorting data
print("\nsorting data...")
Lambda_sort_N4 = np.sort(real_lambda_N4.real)
Lambda_sort_N5 = np.sort(real_lambda_N5.real)
Lambda_sort_N8 = np.sort(real_lambda_N8.real)

print("\ncalculating correlation...")
dif4=dif5=dif8 = 0
for i in range(len(Lambda_sort_N4)-1):
     dif4 = np.append(dif4, Lambda_sort_N4[i+1]-Lambda_sort_N4[i])

for i in range(len(Lambda_sort_N5)-1):
     dif5 = np.append(dif5, Lambda_sort_N5[i+1]-Lambda_sort_N5[i])

for i in range(len(Lambda_sort_N8)-1):
     dif8 = np.append(dif8, Lambda_sort_N8[i+1]-Lambda_sort_N8[i])



#saving to csv
print("\nsaving to csv...")
df_sort4 = pd.DataFrame(Lambda_sort_N4)
df_sort5 = pd.DataFrame(Lambda_sort_N5)
df_sort8 = pd.DataFrame(Lambda_sort_N8)

df_dif4 = pd.DataFrame(dif4)
df_dif5 = pd.DataFrame(dif5)
df_dif8 = pd.DataFrame(dif8)

df_sort4.to_csv("Data/sorted_N4.txt")
df_sort5.to_csv("Data/sorted_N5.txt")
df_sort8.to_csv("Data/sorted_N8.txt")

df_dif4.to_csv("Data/dif4.txt")
df_dif5.to_csv("Data/dif5.txt")
df_dif8.to_csv("Data/dif8.txt")
