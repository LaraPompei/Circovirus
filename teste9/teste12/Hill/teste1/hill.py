import numpy as np
import matplotlib
import matplotlib.pyplot as plt

k1 = 1.0
k2 = 6 
k3 = 20.0

x = np.linspace(0,100,500)

y = k1* (x**k2 / (k3**k2 + x**k2))

plt.figure(figsize=(8,5))
plt.plot(x,y, label=r'$y = k_1\left(\frac{x^{k_2}}{k_3^{k_2} + x^{k_2}}\right)$')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Gráfico da Equação de Hill')
plt.grid(True)
plt.legend()
plt.show()

