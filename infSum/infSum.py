import matplotlib.pyplot as plt
import numpy as np

DTYPE = np.float64
PI = np.pi

T2, T1 = 100.0, 200.0
L, W = 1, 1

N = 20
fact = np.zeros(N)

x, y = 0.5, 0.5
for n in range(1, N):
    fact[n] = fact[n-1] +  (pow(-1,n+1) + 1)/n * np.sin(n*PI*x/L) * np.sinh(n*PI*y/L) /np.sinh(n*PI*W/L) 

fact = np.loadtxt("infSum.txt")
plt.figure()
idx = 1
# plt.scatter(np.linspace(0, N, N)[idx:], T1 + (T2-T1) * 2/PI * fact[idx:])
plt.scatter(np.linspace(0, N, N)[idx:],  fact[idx-1:])
plt.show()
plt.close()