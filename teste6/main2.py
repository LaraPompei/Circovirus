import numpy as np
from scipy.integrate import solve_ivp
import csv

# Constants
eq = 14
t_store = 1000  # Interval of points being saved

# Parameters
#pi_v = 2.38  # Viral replication rate
c_v1 = 32.63  # Maximum viral clearance rate by innate immune system
c_v2 = 8.1e-1  # Half saturation constant

'''
k_v1 = 5.6e-4  # Virus neutralization rate per unit of IgM antibodies
k_v2 = 3.5e-6  # Virus elimination rate per unit of CD8+ T cells
k_v3 = 9.82e-8  # Virus neutralization rate per unit of IgG antibodies
k_v4 = 1.2e-10  # Virus elimination rate per unit of innate immune cells
'''
'''
ajuste 1
pi_v = 1.18
k_v1 = 3.33e-4  # Virus neutralization rate per unit of IgM antibodies
k_v2 = 1.33e-6  # Virus elimination rate per unit of CD8+ T cells
k_v3 = 3.92e-8  # Virus neutralization rate per unit of IgG antibodies
k_v4 = 3.41e-11  # Virus elimination rate per unit of innate immune cells
'''
'''
ajuste 1
pi_v = 2.06

k_v1 = 1.244e-4  # Virus neutralization rate per unit of IgM antibodies
k_v2 = 2.051e-6  # Virus elimination rate per unit of CD8+ T cells
k_v3 = 4.375e-8  # Virus neutralization rate per unit of IgG antibodies
k_v4 = 1.147e-10  # Virus elimination rate per unit of innate immune cells
'''
'''
ajuste 1
pi_v = 2.58

k_v1 = 4.41e-04  # Virus neutralization rate per unit of IgM antibodies
k_v2 = 1.56e-06  # Virus elimination rate per unit of CD8+ T cells
k_v3 = 2.37e-08  # Virus neutralization rate per unit of IgG antibodies
k_v4 = 3.5e-11  # Virus elimination rate per unit of innate immune cells
'''
'''
pi_v = 2.35
k_v1 = 8.5e-5  # Virus neutralization rate per unit of IgM antibodies
k_v2 = 2.32e-6  # Virus elimination rate per unit of CD8+ T cells
k_v3 = 6.16e-8  # Virus neutralization rate per unit of IgG antibodies
k_v4 = 1.88e-11  # Virus elimination rate per unit of innate immune cells
'''

pi_v = 2.275
k_v1 = 3.5e-5  # Virus neutralization rate per unit of IgM antibodies
k_v2 = 2.32e-6  # Virus elimination rate per unit of CD8+ T cells
k_v3 = 6.16e-8  # Virus neutralization rate per unit of IgG antibodies
k_v4 = 1.88e-11  # Virus elimination rate per unit of innate immune cells

#c_pl1 = 1.69e-1
#c_ps1 = 1.44e-2
#delta_IgM = 2.049e-1
#delta_IgG = 1.26e-1


alpha_l = 2.3  # Homeostasis rate of innate immune cells
beta_l = 5.2e-5  # Decay rate of innate immune cells due to virus encounter
delta_l = 1.6e-4  # Natural decay rate of innate immune cells

alpha_ap = 1.5e-3  # Homeostasis rate of immature APCs
beta_ap = 1e-1  # Maturation rate of APCs
c_ap1 = 12e-1  # Maximum maturation rate of APCs
c_ap2 = 2e2  # Half activation constant

delta_apm = 3.1e-1  # Decay rate of mature APCs
alpha_th = 2.17e-4  # Homeostasis rate of CD4+ T cells
beta_th = 1e-7  # Replication rate of naive CD4+ T cells
pi_th = 1e-8  # Replication rate of effector CD4+ T cells
delta_th = 2.2e-1  # Decay rate of effector CD4+ T cells

alpha_tk = 2.17e-4  # Homeostasis rate of CD8+ T cells
beta_tk = 1e-5  # Activation rate of naive CD8+ T cells
pi_tk = 1e-8  # Replication rate of effector CD8+ T cells
delta_tk = 3e-4  # Decay rate of effector CD8+ T cells

alpha_b = 5.5  # Homeostasis rate of B cells
pi_b1 = 4.83e-7  # Activation rate of T-independent B cells
pi_b2 = 1.27e-6  # Activation rate of T-dependent B cells

beta_ps = 1.5e-4  # Differentiation rate of active B cells into short-lived plasma cells
beta_pl = 6.61e-3  # Differentiation rate of active B cells into long-lived plasma cells
beta_bm = 1e-6  # Differentiation rate of active B cells into memory B cells

delta_ps = 2.61e-1  # Decay rate of short-lived plasma cells
delta_pl = 3.2e-3  # Decay rate of long-lived plasma cells
gama_bm = 4.75e-5  # Differentiation rate of memory B cells into long-lived plasma cells

pi_bm1 = 1e-4  # Proliferation rate of memory B cells
pi_bm2 = 2.5e3  # Maximum growth constant

#pi_ps = 19e-4  # Antibody secretion rate per unit of short-lived plasma cells

c_ps1 = 2.03e-2
c_ps2 = 4.5      #dia da subida
c_ps3 = 40
c_ps4 = 40.0
c_ps5 = 5.0

#pi_pl = 1.0e-5  # Antibody secretion rate per unit of long-lived plasma cells

c_pl1 = 6.52732478e-4
c_pl2 = 8.5    
c_pl3 = 60.0
c_pl4 = 80.0
c_pl5 = 21.0

#delta_IgM = 4.42e-1  # Decay rate of IgM
delta_IgM = 5.42e-1  # Decay rate of IgM
delta_IgG = 3.39e-1  # Decay rate of IgG

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

IGG = 1000#54173.4059
IGM = 1000#57361.623 

# Differential equations
def f(t,y):
    dydt = np.zeros(eq)
	#Viremia
    dydt[0] = pi_v * y[0] - ((k_v1 * y[0] * y[11])/(IGM+y[11]))- k_v2 * y[0] * y[6] - ((k_v3 * y[0] * y[12])/(IGG+y[12])) - k_v4 * y[0] * y[13]
    #Apresentadora não madura
    dydt[1] = alpha_ap * (Ap0 - y[1]) - beta_ap * y[1] * (c_ap1 * y[0] / (c_ap2 + y[0]))
    #Apresentadora madura
    dydt[2] = beta_ap * y[1] * (c_ap1 * y[0] / (c_ap2 + y[0])) - delta_apm * y[2]
    #T helper naive
    dydt[3] = alpha_th * (Thn0 - y[3]) - beta_th * y[2] * y[3]
    #T helper efetora
    dydt[4] = beta_th * y[2] * y[3] + pi_th * y[2] * y[4] - delta_th * y[4]
    #T killer naive
    dydt[5] = alpha_tk * (Tkn0 - y[5]) - beta_tk * y[2] * y[5]
    #T killer efetora
    dydt[6] = beta_tk * y[2] * y[5] + pi_tk * y[2] * y[6] - delta_tk * y[6]
    #Celula B
    dydt[7] = alpha_b * (B0 - y[7]) + pi_b2 * y[4] * y[7] - beta_ps * y[2] * y[7] - beta_pl * y[4] * y[7] - beta_bm * y[4] * y[7]
    #Plasma short
    dydt[8] = beta_ps * y[2] * y[7] - delta_ps * y[8]
    #Plasma long
    dydt[9] = beta_pl * y[4] * y[7] - delta_pl * y[9] + gama_bm * y[10]
    #Célula B de memória
    dydt[10] = beta_bm * y[4] * y[7] + pi_bm1 * y[10] * (1 - y[10] / pi_bm2) - gama_bm * y[10]
    #IgM
    dydt[11] = c_ps1 * y[8]  - delta_IgM * y[11]
    #IgG
    dydt[12] = c_pl1 * y[9] - delta_IgG * y[12]
    
	#Tentativa 1 de criar um sistema imune inato genérico
    dydt[13] = alpha_l * (L0 - y[13]) - beta_l * y[13] * y[0] - delta_l * y[13]
    return dydt

# Initial Conditions Vector
y0 = np.array([V0, Ap0, Apm0, Thn0, The0, Tkn0, Tke0, B0, Ps0, Pl0, Bm0, A0, A0, L0])

# Time range
t0 = 0
t_final = 100
h = 0.001
t_eval = np.linspace(t0, t_final, int((t_final-t0)/h))  

# Solve the system using solve_ivp
solution = solve_ivp(f, [t0, t_final], y0, method='LSODA', t_eval=t_eval)

# Extract the results
t_values = solution.t
y_values = solution.y.T  # Transpose to match expected shape

# Generate the CSV file
filename = "output_python.csv"
labels = ['Time', 'V', 'Ap', 'Apm', 'Thn', 'The', 'Tkn', 'Tke', 'B', 'Ps', 'Pl', 'Bm', 'IgM', 'IgG', 'L']

with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(labels)  # Write header
    for i in range(len(t_values)):
        writer.writerow([t_values[i]] + list(y_values[i]))  # Write time and populations row by row

print(f"Data successfully saved to {filename}")
