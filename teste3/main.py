import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.integrate import ode

####################################################################################
# Parâmetros
#####################################################################################
# pi_v = 1.8e-1      # Taxa de replicacao viral
pi_v = 1.61         # Taxa de replicacao viral

# c_v1 = 2.63       # Taxa de clareamento viral maximo pelo sistema inato
c_v1 = 32.63        # Taxa de clareamento viral maximo pelo sistema inato
# c_v2 = 6e-1       # Constante de meia saturacao
c_v2 = 8.1e-1       # Constante de meia saturacao

# k_v1 = 4.82e-5    # Taxa de neutralizacao do virus por unidade anticorpos neutralizantes
k_v1 = 1.6e-4         # Taxa de neutralizacao do virus por unidade anticorpos neutralizantes
# k_v2 = 7.48e-7    # Taxa de eliminacao do virus por unidade de celulas T CD8+
k_v2 = 4e-7     # Taxa de eliminacao do virus por unidade de celulas T CD8+
k_v3 = 4.82e-7      # Taxa de eliminacao do virus por unidade do sistema imune inato
k_v4 = 1.2e-7

alpha_l = 2.3       # taxa de homeostase das celulas do sistema imune inato
beta_l = 2.1e-1    # taxa de decaimento das celulas do sistema imune inato por encontro com virus
delta_l = 1.6e-1   # taxa de decaimento natural das celulas do sistema imune inato

# alpha_ap = 2.5e-3  # Taxa de homeostase das APCs imaturas
alpha_ap = 1.5e-3     # Taxa de homeostase das APCs imaturas
# beta_ap = 5.5e-1    # Taxa de maturacao das APCs 
beta_ap = 1e-1        # Taxa de maturacao das APCs

# c_ap1 = 8e-1        # Taxa de maturacao maxima das APCs
c_ap1 = 12e-1           # bem sensivel - desloca a curva dos anticorpos no eixo x
# c_ap2 = 4e1         # Constante de meia ativacao
c_ap2 = 2e2           # mexe no pico dos anticorpos

# delta_apm = 5.38e-1 # Taxa de morte das APCs maduras
delta_apm = 3.1e-1    # Taxa de morte das APCs maduras
# alpha_th = 2.17e-4  # Taxa de gineistase das celulas T CD4+
alpha_th = 2.17e-4    # Taxa de gineistase das celulas T CD4+
# beta_th = 1e-7      # Taxa de replicacao das celulas T CD4+ naive
beta_th = 1e-7        # Taxa de replicacao das celulas T CD4+ naive
# pi_th = 1e-8        # Taxa de replicacao das celulas T CD4+ efetoras
pi_th = 1e-8          # Taxa de replicacao das celulas T CD4+ efetoras
# delta_th = 2.2e-1   # Taxa de morte das celulas T CD4+ efetoras
delta_th = 2.2e-1     # Taxa de morte das celulas T CD4+ efetoras

# alpha_tk = 2.17e-4  # Taxa de homeostase das celulas T CD8+
alpha_tk = 2.17e-4    # Taxa de homeostase das celulas T CD8+
beta_tk = 1e-5      # Taxa de ativacao das celulas T CD8+ naive
# pi_tk = 1e-8        # Taxa de replicacao das celulas T CD8+ efetoras
pi_tk = 1e-8          # Taxa de replicacao das celulas T CD8+ efetoras
# delta_tk = 3e-4     # Taxa de morte das celulas T CD8+ efetoras
delta_tk = 3e-4       # Taxa de morte das celulas T CD8+ efetoras

# alpha_b = 6.0       # Taxa de homeostase das celulas B
alpha_b = 5.5         # Taxa de homeostase das celulas B

# pi_b1 = 4.83e-6     # Taxa de ativacao das celulas B T-independente
pi_b1 = 4.83e-7       # Taxa de ativacao das celulas B T-independente
# pi_b2 = 1.27e-8     # Taxa de ativacao das celulas B T-dependentes
pi_b2 = 1.27e-8       # Taxa de ativacao das celulas B T-dependentes

# beta_ps = 6.72e-4   # Taxa de diferenciacao das celulas B ativas em plasmocitos de vida curta
beta_ps = 1.85e-2    # Taxa de diferenciacao das celulas B ativas em plasmocitos de vida curta
# beta_pl = 5.61e-6   # Taxa de diferenciacao das celulas B ativas em plasmocitos de vida longa
beta_pl = 6.61e-3     # Taxa de diferenciacao das celulas B ativas em plasmocitos de vida longa

# beta_bm = 1e-6      # Taxa de diferenciacao das celulas B ativas em celulas B de memoria
beta_bm = 1e-6        # afeta a velocidade que o virus decai
# delta_ps = 2.0      # Taxa de morte dos plasmocitos de vida curta
delta_ps = 2.01       # Taxa de morte dos plasmocitos de vida curta
# delta_pl = 2.4e-4   # Taxa de morte dos plasmocitos de vida longa
delta_pl = 3.2e-3     # Taxa de morte dos plasmocitos de vida longa
# gama_bm = 9.75e-4   # Taxa de diferenciacao das celulas B de memoria em plasmocitos de vida longa
gama_bm = 4.75e-1     # Taxa de diferenciacao das celulas B de memoria em plasmocitos de vida longa

pi_bm1 = 1e-5         # Taxa de proliferacao das celulas B de memoria
pi_bm2 = 2.5e3        # Constante de crescimento maximo
# pi_ps = 2e-3        # Taxa de secrecao de anticorpos por unidade de plasmocitos de vida curta

pi_ps = 3.5e-4        # Taxa de secrecao de anticorpos por unidade de plasmocitos de vida curta

c_ps1 = 2.33e-2
c_ps2 = 23.2
c_ps3 = 50.0
c_ps4 = 50.0
c_ps5 = 5.0

c_pl1 = 2.85e-4
c_pl2 = 44.5
c_pl3 = 60.0
c_pl4 = 80.0
c_pl5 = 21.0

# pi_pl = 6.8e-4      # Taxa de secrecao de anticorpos por unidade de plasmocitos de vida longa
pi_pl = 0.0          # Taxa de secrecao de anticorpos por unidade de plasmocitos de vida longa
# delta_a = 4e-2      # Taxa de morte de anticorpos

delta_IgM = 4.42e-1     # Taxa de morte de IgM
delta_IgG = 3.39e-1     # Taxa de morte de IgGV

V0 = 9570.81
L0 = 9570.81
Ap0 = 9.04e6
#Ap0 = 0.83e6
#Ap0 = 10.6e6
Apm0 = 0.0
Thn0 = 1.56e6
#Thn0 = 12.8e6
The0 = 0.0
Tkn0 = 0.91e6
#Tkn0 = 15.4e6  # 28.8e6
Tke0 = 0.0
B0 = 1.39e6
#B0 = 12.4e6
Ps0 = 0.0
Pl0 = 0.0
Bm0 = 0.0
A0 = 0.0

####################################################################################
# Equações
####################################################################################
def V(y, t):
    return pi_v * y[0] - k_v1 * y[0] * y[11] - k_v2 * y[0] * y[6] - k_v3 * y[0] * y[13]

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
    return c_pl1*y[9]*(1 - np.exp(-((t/c_pl2)**c_pl3)) * np.exp(-((t/c_pl4)**c_pl5))) - delta_IgG*y[12]

def L(y, t):
    return alpha_l*(L0 - y[13]) - beta_l*y[13]*(c_ap1*y[0]/(c_ap2 + y[0])) - delta_l*y[13]

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
h = 0.01
num_steps = 1000

# Resolver as EDOs usando o método de Runge-Kutta de 5ª ordem
t_values, y_values = runge_kutta_5(f, t0, y0, h, num_steps)

# Plotar os resultados
labels = ['V', 'Ap', 'Apm', 'Thn', 'The', 'Tkn', 'Tke', 'B', 'Ps', 'Pl', 'Bm', 'IgM', 'IgG', 'L']

for i in range(y_values.shape[1]):
    plt.plot(t_values, y_values[:, i], label=labels[i])

plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()
plt.show()
