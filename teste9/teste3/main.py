import matplotlib.pyplot as plt
from scipy import integrate
from scipy.integrate import ode
import numpy as np
import math

# Constants
eq = 14
t_store = 1000  # Interval of points being saved

# Parameters
pi_v = 3.3  # Viral replication rate
c_v1 = 32.63  # Maximum viral clearance rate by innate immune system
c_v2 = 8.1e-1  # Half saturation constant

k_v1 = 5.6e-4  # Virus neutralization rate per unit of IgM antibodies
k_v2 = 3.5e-6  # Virus elimination rate per unit of CD8+ T cells
k_v3 = 9.82e-8  # Virus neutralization rate per unit of IgG antibodies
k_v4 = 1.2e-10  # Virus elimination rate per unit of innate immune cells

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

beta_ps = 1.55e-2  # Differentiation rate of active B cells into short-lived plasma cells
beta_pl = 6.61e-3  # Differentiation rate of active B cells into long-lived plasma cells
beta_bm = 1e-6  # Differentiation rate of active B cells into memory B cells

delta_ps = 1.81  # Decay rate of short-lived plasma cells
delta_pl = 3.2e-3  # Decay rate of long-lived plasma cells
gama_bm = 4.75e-5  # Differentiation rate of memory B cells into long-lived plasma cells

pi_bm1 = 1e-4  # Proliferation rate of memory B cells
pi_bm2 = 2.5e3  # Maximum growth constant

pi_ps = 19e-4  # Antibody secretion rate per unit of short-lived plasma cells

c_ps1 = 1.63e-2
c_ps2 = 22
c_ps3 = 50.0
c_ps4 = 40.0
c_ps5 = 5.0

pi_pl = 1.0e-5  # Antibody secretion rate per unit of long-lived plasma cells

c_pl1 = 2.85e-4
c_pl2 = 44.5
c_pl3 = 60.0
c_pl4 = 80.0
c_pl5 = 21.0

delta_IgM = 4.42e-1  # Decay rate of IgM
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

####################################################################################
# Equações
####################################################################################
def V(y, t):
    return pi_v * y[0] - k_v1 * y[0] * y[11] - k_v2 * y[0] * y[6] - k_v3 * y[0] * y[12] - k_v4*y[0]*y[13]

def Ap(y, t):
    return alpha_ap * (Ap0 - y[1]) - beta_ap * y[1] * (c_ap1 * y[0] / (c_ap2 + y[0]))

def Apm(y, t):
    return beta_ap * y[1] * (c_ap1 * y[0] / (c_ap2 + y[0])) - delta_apm * y[2]

def Thn(y, t):
    return alpha_th * (Thn0 - y[3]) - beta_th * y[2] * y[3]

def The(y, t):
    return beta_th * y[2] * y[3] + pi_th * y[2] * y[4] - delta_th * y[4]

def Tkn(y, t):
    return alpha_tk * (Tkn0 - y[5]) - beta_tk * y[2] * y[5]

def Tke(y, t):
    return beta_tk * y[2] * y[5] + pi_tk * y[2] * y[6] - delta_tk * y[6]

def B(y, t):
    return alpha_b * (B0 - y[7]) + pi_b2 * y[4] * y[7] - beta_ps * y[2] * y[7] - beta_pl * y[4] * y[7] - beta_bm * y[4] * y[7]

def Ps(y, t):
    return beta_ps * y[2] * y[7] - delta_ps * y[8]

def Pl(y, t):
    return beta_pl * y[4] * y[7] - delta_pl * y[9] + gama_bm * y[10]

def Bm(y, t):
    return beta_bm * y[4] * y[7] + pi_bm1 * y[10] * (1 - y[10] / pi_bm2) - gama_bm * y[10]
def IgM(y, t):
    return c_ps1*y[8]*(1 - np.exp(-((t/c_ps2)**c_ps3)) * np.exp(-((t/c_ps4)**c_ps5))) - delta_IgM*y[11]

def IgG(y, t):
    return c_pl1*y[9] * (1 - np.exp(-((t/c_pl2)**c_pl3)) * np.exp(-((t/c_pl4)**c_pl5))) - delta_IgG*y[12]

def L(y, t):
    return alpha_l*(L0 - y[13]) - beta_l*y[13]*y[0] - delta_l*y[13]

# Resolvendo o sistema de EDOs
def f(t, y):
    return np.array([V(y, t), Ap(y, t), Apm(y, t), Thn(y, t), The(y, t), Tkn(y, t), Tke(y, t), B(y, t), Ps(y, t), Pl(y, t), Bm(y, t), IgM(y, t), IgG(y, t), L(y, t)])

def runge_kutta_5(f, t0, y0, h, num_steps):
    t_values = [t0]
    y_values = [y0]

    for _ in range(num_steps):
        t = t_values[-1]
        y = y_values[-1]

        k1 = h * f(t, y)
        k2 = h * f(t + h / 4, y + k1 / 4)
        k3 = h * f(t + 3 * h / 8, y + 3 * k1 / 32 + 9 * k2 / 32)
        k4 = h * f(t + 12 * h / 13, y + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197)
        k5 = h * f(t + h, y + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104)
        k6 = h * f(t + h / 2, y - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40)

        t_values.append(t + h)
        y_values.append(y + 25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - k5 / 5)

    return np.array(t_values), np.array(y_values)

# Exemplo de uso

# Definir condições iniciais e parâmetros
t0 = 0
y0 = np.array([V0, Ap0, Apm0, Thn0, The0, Tkn0, Tke0, B0, Ps0, Pl0, Bm0, A0, A0, L0])
h = 0.001
t_final = 150

num_steps = t_final*1000

# Solve the system using the improved Runge-Kutta method
t_values, y_values = runge_kutta_5(f, t0, y0, h, num_steps)

# Plotar os resultados
labels = ['V', 'Ap', 'Apm', 'Thn', 'The', 'Tkn', 'Tke', 'B', 'Ps', 'Pl', 'Bm', 'IgM', 'IgG', 'L']
# Plotting each population separately
fig, axs = plt.subplots(len(labels), figsize=(10, 30))
fig.tight_layout(pad=4.0)

for i in range(len(labels)):
    plt.figure()  # Create a new figure
    plt.plot(t_values, y_values[:, i], label=labels[i])
    plt.xlabel('Time (days)')
    plt.ylabel(labels[i])
    plt.legend()
    plt.title(f'Time evolution of {labels[i]}')
    plt.show()  # Display the plot
