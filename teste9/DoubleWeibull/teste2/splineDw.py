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

k1 = 2.0e6
k2 = 2e6 #*0.01* x.max()
k3 = 80.0
k4 = 6e6 #*0.1*x.max()
k5 = 40e6 #sem muito efeito 

y = k1 * (1 - np.exp(- (x / k2)**k3)) * np.exp(- (x / k4)**k5)
 
plt.figure(figsize=(8,5))
#plt.plot(t, x, label='data spline')
plt.plot(t, y, label=r'$y = k_1\left(1 - e^{-(x/k2)^{k3}}\right)e^{-(x/k4)^k5}$')
plt.xlabel('tempo')
plt.ylabel('viremia')
plt.title('Gráfico da Equação de Double Weilbull com Apresentadora')
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
