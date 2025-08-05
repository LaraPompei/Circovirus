import numpy as np
import matplotlib
import matplotlib.pyplot as plt

k1 = 1.0
k2 = 15.0 
k3 = 10.0
k4 = 50.0 
k5 = 5.0

x = np.linspace(0,100,500)

y = k1 * (1 - np.exp(- (x / k2)**k3)) * np.exp(- (x / k4)**k5)

plt.figure(figsize=(8,5))
plt.plot(x,y, label=r'$y = k_1\left(1 - e^{-(x/k2)^{k3}}\right)e^{-(x/k4)^k5}$')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Gráfico da Equação de DW')
plt.grid(True)
plt.legend()
plt.show()

