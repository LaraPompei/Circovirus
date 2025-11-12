import numpy as np
from scipy.integrate import solve_ivp
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt

#dados 
target_antibody_data = np.genfromtxt('./dados/plot-IGM.csv', delimiter=',',skip_header=1)

dados_viremia = pd.read_csv('./dados/viremiaPorcoInoculado.csv')

#PARAMETROS
n = 70 
eq = 13+n
t_store = 1000

alpha_l = 2.3  # Homeostasis rate of innate immune cells
beta_l = 5.2e-5  # Decay rate of innate immune cells due to virus encounter
delta_l = 1.6e-4  # Natural decay rate of innate immune cells

#alpha_ap = 1.5e-3  # Homeostasis rate of immature APCs
alpha_ap = 1.5e-2
#beta_ap = 1e-3  # Maturation rate of APCs
beta_ap = 9e-1
tau = 21 

c_ap1 = 1       # Maximum maturation rate of APCs
c_ap2 = 40      # Half activation constant
c_ap3 = 10e6     # Day the APm start being produced

delta_api = 6.1e-2  # Decay rate of intermediary APC

#gama_api = 2*3.1e-5   #Maturation rate of intermediary APCs
gama_api = 2*3.1e-1

delta_apm = 5.1e-1  #Decay rate of mature APCs 
#taum = 20

#alpha_th = 2.17e-4  # Homeostasis rate of CD4+ T cells
alpha_th = 7.17e-2

#beta_th = 1e-5  # Replication rate of naive CD4+ T cells
beta_th = 1e-6
pi_th = 1e-7  # Replication rate of effector CD4+ T cells
delta_th = 2.2e0  # Decay rate of effector CD4+ T cells

#alpha_tk = 2.17e-4  # Homeostasis rate of CD8+ T cells
alpha_tk = 2.17e-1  # Homeostasis rate of CD8+ T cells

#beta_tk = 1e-7  # Activation rate of naive CD8+ T cells
beta_tk = 1e-8
#pi_tk = 1e-8
pi_tk = 1e-8  # Replication rate of effector CD8+ T cells
#delta_tk = 3e-4  # Decay rate of effector CD8+ T cells
delta_tk = 3e-1

#alpha_b = 5.5e1  # Homeostasis rate of B cells
alpha_b = 5.5e1
pi_b1 = 4.83e-7  # Activation rate of T-independent B cells
pi_b2 = 1.27e-6  # Activation rate of T-dependent B cells


#beta_ps = 4.5e-4  # Differentiation rate of active B cells into short-lived plasma cells
beta_ps = 4.5e-5

#beta_pl = 6.61e-3  # Differentiation rate of active B cells into long-lived plasma cells
beta_pl = 0.0#6.61e-12  # Differentiation rate of active B cells into long-lived plasma cells

beta_bm = 0#1e-6  # Differentiation rate of active B cells into memory B cells

delta_ps = 7.61e1  # Decay rate of short-lived plasma cells

delta_pl = 3.2e-10  # Decay rate of long-lived plasma cells
gama_bm = 4.75e-5  # Differentiation rate of memory B cells into long-lived plasma cells

pi_bm1 = 1e-4  # Proliferation rate of memory B cells
pi_bm2 = 2.5e3  # Maximum growth constant

#pi_ps = 19e-4  # Antibody secretion rate per unit of short-lived plasma cells
#c_ps1 = 2.03e-2
#c_ps1 = 5.8e7
#c_ps2 = 15      #dia da subida
c_ps1 = 1.24e2#33635.7168859# 1.51e6 
c_ps2 = 2.3e-1
c_ps3 = 9.03e4

#c_ps4 = 40.0
#c_ps5 = 5.0

c_pl1 = 5.6e-25
c_pl2 = 50.5    #dia da subida
c_pl3 = 60.0
c_pl4 = 80.0
c_pl5 = 21.0

#delta_IgM = 4.42e-1  # Decay rate of IgM
delta_IgM = 317.89443611#5.42e1  # Decay rate of IgM
#delta_IgG = 3.39e-1  # Decay rate of IgG
delta_IgG = 3.39e3  # Decay rate of IgG

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



def model(params, t):
    delta_apm, c_ps1, delta_IgM, beta_ps, pi_b2, alpha_b, beta_th, pi_th, alpha_th = params
    #print(delta_apm, c_ps1, delta_IgM, beta_ps, pi_b2, alpha_b, beta_th, pi_th, alpha_th)
    y0 = np.array([V0, Ap0, Apm0, Thn0, The0, Tkn0, Tke0, B0, Ps0, Pl0, Bm0, A0, A0])
    for k in range(eq-n, eq):
        y0 = np.append(y0,Apm0)
    #print(f"len(y0) = {y0.size}")
    dados_viremia.columns = dados_viremia.columns.str.strip()

    told = dados_viremia['x'].values
    xold = dados_viremia['y'].values

    spline = CubicSpline(told, xold)
    spline_x = spline.derivative(1)
    # Differential equations
    def f(t,y):
        dydt = np.zeros(eq)
        #virus
        dydt[0] = 0 if t > 54 else (spline_x(t))
        #Apresentadora naive
        dydt[1] = alpha_ap * (Ap0 - y[1]) - beta_ap * y[1] * ((c_ap1 * y[0]**c_ap2) / (c_ap3**c_ap2 + y[0]**c_ap2))
        #Apresentadora ideal
        dydt[2] = beta_ap * y[1] * ((c_ap1 * y[0]**c_ap2) / (c_ap3**c_ap2 + y[0]**c_ap2)) - delta_apm * y[2]
        #Apresentadora intermediaria 
        dydt[13] = (n*(y[2] - y[13]))/tau
        #Apresentadora Madura
        for i in range(eq-n+1,eq):
            #print(i-1,i)
            dydt[i] = (n * (y[i-1] - y[i]))/tau 
        #T helper naive
        dydt[3] = alpha_th * (Thn0 - y[3]) - beta_th * y[eq-1] * y[3]
        #T helper efetora 
        dydt[4] = beta_th * y[eq-1] * y[3] + pi_th * y[eq-1] * y[4] - delta_th * y[4]
        #T killer naive
        dydt[5] = alpha_tk * (Tkn0 - y[5]) - beta_tk * y[eq-1] * y[5]
        #T killer efetora 
        dydt[6] = beta_tk * y[eq-1] * y[5] + pi_tk * y[eq-1] * y[6] - delta_tk * y[6]
        #B
        dydt[7] = alpha_b * (B0 - y[7]) + pi_b2 * y[4] * y[7] - beta_ps * y[eq-1] * y[7] - beta_pl * y[4] * y[7] - beta_bm * y[4] * y[7]
        #Ps
        dydt[8] = beta_ps * y[eq-1] * y[7] - delta_ps * y[8]
        #Pl
        dydt[9] = beta_pl * y[4] * y[7] - delta_pl * y[9] + gama_bm * y[10]
        #Bm
        dydt[10] = beta_bm * y[4] * y[7] + pi_bm1 * y[10] * (1 - y[10] / pi_bm2) - gama_bm * y[10]
        #IgM
        dydt[11] = c_ps1* y[8] - delta_IgM * y[11]
        #IgG 
        dydt[12] = 0#c_pl1 * y[9] - delta_IgG * y[12]
        if np.any(np.isnan(dydt)) or np.any(np.isinf(dydt)):
            return np.full_like(dydt, -1e6)
        return dydt

    #Solve the system
    solution = solve_ivp(f, [t0, t_final], y0, method='LSODA', t_eval=t)
    return solution.t, solution.y.T

def erro(params):
    print("Testing params:", params)
    #solving the model
    t_values, y_values = model(params, t_eval)

    #define antibody values
    model_antibody = y_values[:,11]
    if np.any(np.isnan(model_antibody)) or np.any(np.isinf(model_antibody)):
            return 1e20
    
    #defining time and IgM population values for the experimental data
    t_exp = target_antibody_data[:,0]
    IgM_exp = target_antibody_data[:,1]

    #interpolate the t and IgM values from the model
    model_antibody_interp = interp1d(t_values, model_antibody, kind = 'linear', fill_value='extrapolate')
    
    #define IgM in the experimental time
    IgM = model_antibody_interp(t_exp) 
    
    l2_error_ant = np.linalg.norm(IgM - IgM_exp)/np.linalg.norm(IgM_exp)

    total_error = l2_error_ant*100

    print(f"total_error = {total_error}")
    
    IgM_model = y_values[:, 1]        # IgM do modelo
    t_exp = target_antibody_data[:, 0] # tempos experimentais
    IgM_exp = target_antibody_data[:, 1] # IgM experimental

    print(f"IgM_model min={IgM_model.min():.2e}, max={IgM_model.max():.2e}")

    plt.figure(figsize=(6,4))
    plt.plot(t_values, IgM_model, label='IgM (modelo)', color='blue')
    #plt.scatter(t_exp, IgM_exp, color='red', label='IgM (experimental)', zorder=3)
    plt.xlabel('Tempo (dias)')
    plt.ylabel('Concentração')
    plt.title('Comparação: IgM do modelo vs dados experimentais')
    plt.legend()
    plt.tight_layout()
    plt.show()


    return total_error

t0 = 0
t_final =100 
h = 0.001
t_eval = np.linspace(t0, t_final, int((t_final-t0)/h))

# Initial guesses and bounds for parameters
n_bounds = (0,50) 
tau_bounds = (0,50)
delta_apm_bounds = (5e-2,5e0)
c_ps1_bounds = (1.0e1,1.0e3)
c_ps2_bounds = (0,50)
c_ps3_bounds =  (0,240000)
delta_IgM_bounds = (3.17e1,1.17e3)
beta_ps_bounds =  (4.5e-6,4.5e-4)
pi_b2_bounds = (1.27e-7,1.27e-5)
alpha_b_bounds = (5.5e0,5.5e2)
beta_th_bounds = (1e-7,1e-5)
pi_th_bounds = (1e-8,1e-6)
alpha_th_bounds = (1e-3,1e-1)

bounds = [
    delta_apm_bounds,
    c_ps1_bounds,  
    delta_IgM_bounds, 
    beta_ps_bounds,  
    pi_b2_bounds, 
    alpha_b_bounds, 
    beta_th_bounds, 
    pi_th_bounds, 
    alpha_th_bounds  
]
result = differential_evolution(erro, bounds, strategy='best1bin', maxiter=1000, popsize=15, init='latinhypercube', tol=0.01, mutation=(0.5, 1), recombination=0.9, disp=True, x0 = [delta_apm,c_ps1,delta_IgM, beta_ps, pi_b2, alpha_b, beta_th, pi_th, alpha_th])
print("Baseline error:", erro([delta_apm,c_ps1,delta_IgM,beta_ps,pi_b2,alpha_b,beta_th,pi_th,alpha_th]))

optimal_params = result.x

final_model_output = model(optimal_params, t_eval)

# Generate the CSV file
filename = "output_DE.csv"
labels = ['Time', 'V', 'Ap', 'Apm', 'Thn', 'The', 'Tkn', 'Tke', 'B', 'Ps', 'Pl', 'Bm', 'IgM', 'IgG', 'L']
for i in range(0, n):
    labels.append('Ap'+str(i))

with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(labels)  # Write header
    for i in range(len(t_eval)):
        writer.writerow([t_eval[i]] + list(final_model_output[1][:, i]))  # Write time and populations row by row

print("Best parameters found:")
print(optimal_params)
print("Best error value:")
print(result.fun)

print(f"Data successfully saved to {filename}")
    
