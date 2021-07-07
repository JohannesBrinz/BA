import numpy as np
from scipy import linalg as lg
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import cmath



#data
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
dif12 = pd.read_csv("Data/dif12.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif12 = dif12["number"].to_numpy()
dif16 = pd.read_csv("Data/dif16.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif16 = dif16["number"].to_numpy()
dif20 = pd.read_csv("Data/dif20.txt", sep = ",", header = 0, \
    names = ["index", "number"])
dif20 = dif20["number"].to_numpy()

#counter
counter = 0
for i in range(len(dif3)):
    if dif3[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=3:\t", counter/len(dif3))

counter = 0
for i in range(len(dif4)):
    if dif4[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=4:\t", counter/len(dif4))

counter = 0
for i in range(len(dif5)):
    if dif5[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=5:\t", counter/len(dif5))

counter = 0
for i in range(len(dif6)):
    if dif6[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=6:\t", counter/len(dif6))

counter = 0
for i in range(len(dif7)):
    if dif7[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=7:\t", counter/len(dif7))

counter = 0
for i in range(len(dif8)):
    if dif8[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=8:\t", counter/len(dif8))

counter = 0
for i in range(len(dif12)):
    if dif12[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=12:\t", counter/len(dif12))


counter = 0
for i in range(len(dif16)):
    if dif16[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=16:\t", counter/len(dif16))

counter = 0
for i in range(len(dif20)):
    if dif20[i] < 1e-12:
        counter += 1
print("The fraction of zero distances is for N=20:\t", counter/len(dif20))
