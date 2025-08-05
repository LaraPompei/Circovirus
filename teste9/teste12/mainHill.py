import numpy as np
from scipy.integrate import solve_ivp
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline

#PARAMETROS
eq = 14
t_store = 1000

alpha_l = 2.3  # Homeostasis rate of innate immune cells
beta_l = 5.2e-5  # Decay rate of innate immune cells due to virus encounter
delta_l = 1.6e-4  # Natural decay rate of innate immune cells

alpha_ap = 1.5e-3  # Homeostasis rate of immature APCs
#beta_ap = 1e-3  # Maturation rate of APCs
beta_ap = 1e-1#8e3
taui = 2e2

c_ap1 = 1       # Maximum maturation rate of APCs
c_ap2 = 20      # Half activation constant
c_ap3 = 20     # Day the APm start being produced

delta_api = 3.1e-1  # Decay rate of intermediary APC

#gama_api = 2*3.1e-5   #Maturation rate of intermediary APCs
gama_api = 2e-1

delta_apm = 3.1e-6  #Decay rate of mature APCs 
taum = 2e1

alpha_th = 2.17e-4  # Homeostasis rate of CD4+ T cells
#beta_th = 1e-5  # Replication rate of naive CD4+ T cells
beta_th = 1e1
pi_th = 1e-8  # Replication rate of effector CD4+ T cells
delta_th = 2.2e-1  # Decay rate of effector CD4+ T cells

alpha_tk = 2.17e-4  # Homeostasis rate of CD8+ T cells
#beta_tk = 1e-7  # Activation rate of naive CD8+ T cells
beta_tk = 1e-4
pi_tk = 1e-8  # Replication rate of effector CD8+ T cells
delta_tk = 3e-4  # Decay rate of effector CD8+ T cells

alpha_b = 5.5  # Homeostasis rate of B cells
pi_b1 = 4.83e-7  # Activation rate of T-independent B cells
pi_b2 = 1.27e-6  # Activation rate of T-dependent B cells

#beta_ps = 4.5e-4  # Differentiation rate of active B cells into short-lived plasma cells
beta_ps = 4.5e-3

beta_pl = 6.61e-3  # Differentiation rate of active B cells into long-lived plasma cells
beta_bm = 1e-6  # Differentiation rate of active B cells into memory B cells

delta_ps = 2.61e-1  # Decay rate of short-lived plasma cells
delta_pl = 3.2e-3  # Decay rate of long-lived plasma cells
gama_bm = 4.75e-5  # Differentiation rate of memory B cells into long-lived plasma cells

pi_bm1 = 1e-4  # Proliferation rate of memory B cells
pi_bm2 = 2.5e3  # Maximum growth constant

#pi_ps = 19e-4  # Antibody secretion rate per unit of short-lived plasma cells
#c_ps1 = 2.03e-2
c_ps1 = 2.03e-1
c_ps2 = 23.5      #dia da subida
c_ps3 = 40
c_ps4 = 40.0
c_ps5 = 5.0

c_pl1 = 5.6e-9
c_pl2 = 50.5    #dia da subida
c_pl3 = 60.0
c_pl4 = 80.0
c_pl5 = 21.0

#delta_IgM = 4.42e-1  # Decay rate of IgM
delta_IgM = 5.42e-2  # Decay rate of IgM
delta_IgG = 3.39e-1  # Decay rate of IgG

c_apm1 = 1
c_apm2 = 20
c_apm3 = 2e4

# Initial Conditions
V0 = 9570.81
L0 = 9.04e6
Ap0 = 0.83e6
Apm0 = 0.0
Thn0 = 1.56e6
The0 = 0.0
Tkn0 = 0.91e6
Tke0 = 0.0
B0 = 1.39e6
Ps0 = 0.0
Pl0 = 0.0
Bm0 = 0.0
A0 = 0.0

# Differential equations
def f(t,y,spline_x):
    dydt = np.zeros(eq)
    #virus
    dydt[0] = 0 if t > 54 else (spline_x(t))
    #Apresentadora naive
    dydt[1] = alpha_ap * (Ap0 - y[1]) - beta_ap * y[1] * ((c_ap1 * y[0]**c_ap2) / (c_ap3**c_ap2 + y[0]**c_ap2))
    #Apresentadora intermediaria    
    dydt[2] = (beta_ap * y[1] * ((c_ap1 * y[0]**c_ap2) / (c_ap3**c_ap2 + y[0]**c_ap2)) - y[2])/taui
    #dydt[2] = (beta_ap * y[1] * ((c_ap1 * y[0]**c_ap2) / (c_ap3**c_ap2 + y[0]**c_ap2)) - ((c_apm1 * y[2]**c_apm2) / (c_apm3**c_apm2 + y[2]**c_apm2)))/taui
    #dydt[2] = beta_ap * y[1] * (c_ap1 * y[0]**c_ap2 / (c_ap3**c_ap2 + y[0]**c_ap2)) - gama_api * ((c_apm1 * y[2]**c_apm2) / (c_apm3**c_apm2 + y[2]**c_apm2))
   #Apresentadora Madura 
    dydt[3] = (gama_api *y[2] - y[3])/taum
    #T helper naive
    #dydt[3] = (gama_api* ((c_apm1 * y[2]**c_apm2) / (c_apm3**c_apm2 + y[2]**c_apm2)) - delta_apm * y[3])
    dydt[4] = alpha_th * (Thn0 - y[4]) - beta_th * y[3] * y[4]
    #T helper efetora 
    dydt[5] = beta_th * y[3] * y[4] + pi_th * y[3] * y[5] - delta_th * y[5]
    #T killer naive
    dydt[6] = alpha_tk * (Tkn0 - y[6]) - beta_tk * y[3] * y[6]
    #T killer efetora 
    dydt[7] = beta_tk * y[3] * y[6] + pi_tk * y[3] * y[7] - delta_tk * y[7]
    #dydt[7] = alpha_b * (B0 - y[7]) + pi_b2 * y[4] * y[7] - beta_ps * y[2] * (c_ap1 * y[7]**c_ap2 / (c_ap3**c_ap2 + y[7]**c_ap2)) - beta_pl * y[4] * y[7] - beta_bm * y[4] * y[7]
    #B
    #dydt[8] = beta_ps * y[2] * (c_ap1 * y[7]**c_ap2 / (c_ap3**c_ap2 + y[7]**c_ap2)) - delta_ps * y[8]
    dydt[8] = alpha_b * (B0 - y[8]) + pi_b2 * y[5] * y[8] - beta_ps * y[3] * y[8] - beta_pl * y[5] * y[8] - beta_bm * y[5] * y[8]
    #Ps
    dydt[9] = beta_ps * y[3] * y[8] - delta_ps * y[9]
    #Pl
    dydt[10] = beta_pl * y[5] * y[8] - delta_pl * y[10] + gama_bm * y[11]
    #Bm
    dydt[11] = beta_bm * y[5] * y[8] + pi_bm1 * y[11] * (1 - y[11] / pi_bm2) - gama_bm * y[11]
    #IgM
    dydt[12] = c_ps1 * y[9] - delta_IgM * y[12]
    #IgG 
    dydt[13] = c_pl1 * y[10] - delta_IgG * y[13]
    return dydt

arquivo_viremia = '/home/larapompei/Documents/Disciplinas/Mestrado/Pesquisa/Circovirus/teste7/dados/viremiaPorcoInoculado.csv'
dados_viremia = pd.read_csv(arquivo_viremia)
dados_viremia.columns = dados_viremia.columns.str.strip()

told = dados_viremia['x'].values
xold = dados_viremia['y'].values

spline_x = CubicSpline(told, xold)
#t = np.linspace(told.min(), told.max(), 500)

# Initial Conditions Vector
y0 = np.array([V0, Ap0, Apm0, Apm0, Thn0, The0, Tkn0, Tke0, B0, Ps0, Pl0, Bm0, A0, A0])

# Time range
t0 = 0
t_final =100 
h = 0.001
t_eval = np.linspace(t0, t_final, int((t_final-t0)/h)) 
c_ap3 = spline_x(t_eval).max()*0.8
# Solve the system using solve_ivp
solution = solve_ivp(f, [t0, t_final], y0, method='LSODA', t_eval=t_eval, args=(spline_x.derivative(1),))

t_values = solution.t
y_values = solution.y.T  

filename = "output_python.csv"
labels = ['Time', 'V', 'Ap','Api', 'Apm', 'Thn', 'The', 'Tkn', 'Tke', 'B', 'Ps', 'Pl', 'Bm', 'IgM', 'IgG']

#print(told.max())
with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(labels)  # Write header
    for i in range(len(t_values)):
        writer.writerow([t_values[i]] + list(y_values[i]))  

print(f"Data successfully saved to {filename}")

dydt_values = np.array([f(t, y, spline_x.derivative(1)) for t, y in zip(t_values, y_values)])
dydt0 = dydt_values[:, 0]
dydt1 = dydt_values[:, 1]

#hill = beta_ap * dydt1 * (c_ap1 * dydt0**c_ap2 / (c_ap3**c_ap2 + dydt0**c_ap2))

#plt.figure(figsize=(8, 5))
#plt.plot(t_eval, hill, label='Hill curve')
#plt.grid(True)
#plt.legend()
#plt.title('Hill equation before integration')
#plt.show()


#hill2 = beta_ap * y_values[:,1] * (c_ap1 * y_values[:,0]**c_ap2 / (c_ap3**c_ap2 + y_values[:,0]**c_ap2))

#plt.plot(t_eval, hill2, label='Hill curve')
#plt.grid(True)
#plt.legend()
#plt.title('Hill equation after integration')
#plt.show()
