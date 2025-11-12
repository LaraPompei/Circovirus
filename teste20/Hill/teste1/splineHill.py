import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline

arquivo = '/home/larapompei/Documents/Disciplinas/Mestrado/Pesquisa/Circovirus/teste7/dados/viremiaPorcoInoculado.csv'
dados = pd.read_csv(arquivo)
dados.columns = dados.columns.str.strip()

told = dados['x'].values
xold = dados['y'].values

spline = CubicSpline(told, xold)

t = np.linspace(told.min(), told.max(), 500)
x = spline(t)

k1 = 1.0
k2 = 20 
k3 = 0.8*x.max()

y = k1*x* (x**k2 / (k3**k2 + x**k2))

 
plt.figure(figsize=(8,5))
#plt.plot(t, x, label='data spline')
plt.plot(t, y, label=r'$y = k_1\left(\frac{x^{k_2}}{k_3^{k_2} + x^{k_2}}\right)$')
plt.xlabel('tempo')
plt.ylabel('viremia')
plt.title('Gráfico da Equação de Hill')
plt.grid(True)
plt.legend()
plt.show()

df_result = pd.DataFrame({
    'tempo': t,
    'viremia_spline': x,
    'hill_transform': y
})
output_path = "/home/larapompei/Documents/Disciplinas/Mestrado/Pesquisa/Circovirus/Hill/viremia_hill_result.csv"
df_result.to_csv(output_path, index=False)
output_path
