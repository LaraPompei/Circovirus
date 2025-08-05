import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline

arquivo_viremia = '/home/larapompei/Documents/Disciplinas/Mestrado/Pesquisa/Circovirus/teste7/dados/viremiaPorcoInoculado.csv'
dados_viremia = pd.read_csv(arquivo_viremia)
dados_viremia.columns = dados_viremia.columns.str.strip()

told = dados_viremia['x'].values
xold = dados_viremia['y'].values

spline_x = CubicSpline(told, xold)

t = np.linspace(told.min(), told.max(), 500)
x = spline_x(t)

arquivo_output = 'output.csv'
dados_output = pd.read_csv(arquivo_output)
dados_output.columns = dados_output.columns.str.strip()

t_output = dados_output['t'].values
Apm = dados_output['Apm'].values

spline_z = CubicSpline(t_output, Apm)
z = spline_z(t)

k1 = 2.0
k2 = 20
k3 = x.max()
k3t = 20

y = k1 * z * (x**k2 / (k3**k2 + x**k2))
#hill = k1* (t**k2 / (k3t**k2 + t**k2))*1e5

plt.figure(figsize=(8,5))
#plt.plot(t, x, label='viremia spline(x)')
#plt.plot(t, z, label='Apm spline(z)')
plt.plot(t, y, label=r'$y = k_1 \cdot z \cdot \frac{x^{k_2}}{k_3^{k_2} + x^{k_2}}$', linewidth=2)
#plt.plot(t,hill,label=r'$y = k_1 \cdot \frac{x^{k_2}}{k_3^{k_2} + x^{k_2}} * 10^{5}$')
plt.xlabel('tempo')
plt.ylabel('população')
plt.title('Hill com duas populações diferentes')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

