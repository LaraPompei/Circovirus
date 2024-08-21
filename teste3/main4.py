import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import csv

# Constants
eq = 14
t_store = 1000  # Interval of points being saved

# Parameters (bounds are set later)
k_v1 = 5.6e-4
k_v2 = 3.5e-6
k_v3 = 9.82e-8
k_v4 = 1.2e-10
alpha_l = 2.3
beta_l = 5.2e-5
delta_l = 1.6e-4
c_ap1 = 12e-1
c_ap2 = 2e2
alpha_th = 2.17e-4
beta_th = 1e-7
pi_th = 1e-8
delta_th = 2.2e-1
alpha_tk = 2.17e-4
beta_tk = 1e-5
pi_tk = 1e-8
delta_tk = 3e-4
alpha_b = 5.5
pi_b1 = 4.83e-7
pi_b2 = 1.27e-6
beta_ps = 1.55e-2
beta_pl = 6.61e-3
beta_bm = 1e-6
delta_ps = 1.81
delta_pl = 3.2e-3
gama_bm = 4.75e-5
pi_bm1 = 1e-4
pi_bm2 = 2.5e3
pi_ps = 19e-4
c_ps1 = 1.63e-2
c_ps2 = 22
c_ps3 = 50.0
c_ps4 = 40.0
c_ps5 = 5.0
pi_pl = 1.0e-5
c_pl1 = 2.85e-4
c_pl2 = 44.5
c_pl3 = 60.0
c_pl4 = 80.0
c_pl5 = 21.0
delta_IgM = 4.42e-1
delta_IgG = 3.39e-1

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

# Time definitions
t0 = 0  # Start time
t_final = 150  # End time

# Load target data for antibody population and viremia
target_antibody_data = np.genfromtxt('./dados/anticorposPorcoInoculado.csv', delimiter=',', skip_header=1)
target_viremia_data = np.genfromtxt('./dados/viremiaPorcoInoculado.csv', delimiter=',', skip_header=1)

# Define your model with the parameters to be optimized
def model(params, t_eval):
    pi_v, alpha_ap, beta_ap, delta_apm = params

    def f(t, y):
        dydt = np.zeros(eq)
        dydt[0] = pi_v * y[0] - k_v1 * y[0] * y[11] - k_v2 * y[0] * y[6] - k_v3 * y[0] * y[12] - k_v4 * y[0] * y[13]
        dydt[1] = alpha_ap * (Ap0 - y[1]) - beta_ap * y[1] * (c_ap1 * y[0] / (c_ap2 + y[0]))
        dydt[2] = beta_ap * y[1] * (c_ap1 * y[0] / (c_ap2 + y[0])) - delta_apm * y[2]
        dydt[3] = alpha_th * (Thn0 - y[3]) - beta_th * y[2] * y[3]
        dydt[4] = beta_th * y[2] * y[3] + pi_th * y[2] * y[4] - delta_th * y[4]
        dydt[5] = alpha_tk * (Tkn0 - y[5]) - beta_tk * y[2] * y[5]
        dydt[6] = beta_tk * y[2] * y[5] + pi_tk * y[2] * y[6] - delta_tk * y[6]
        dydt[7] = alpha_b * (B0 - y[7]) + pi_b2 * y[4] * y[7] - beta_ps * y[2] * y[7] - beta_pl * y[4] * y[7] - beta_bm * y[4] * y[7]
        dydt[8] = beta_ps * y[2] * y[7] - delta_ps * y[8]
        dydt[9] = beta_pl * y[4] * y[7] - delta_pl * y[9] + gama_bm * y[10]
        dydt[10] = beta_bm * y[4] * y[7] + pi_bm1 * y[10] * (1 - y[10] / pi_bm2) - gama_bm * y[10]
        dydt[11] = c_ps1 * y[8] * (1 - np.exp(-((t / c_ps2)**c_ps3)) * np.exp(-((t / c_ps4)**c_ps5))) - delta_IgM * y[11]
        dydt[12] = c_pl1 * y[9] * (1 - np.exp(-((t / c_pl2)**c_pl3)) * np.exp(-((t / c_pl4)**c_pl5))) - delta_IgG * y[12]
        dydt[13] = alpha_l * (L0 - y[13]) - beta_l * y[13] * y[0] - delta_l * y[13]
        return dydt

    # Initial Conditions
    y0 = np.array([V0, Ap0, Apm0, Thn0, The0, Tkn0, Tke0, B0, Ps0, Pl0, Bm0, A0, A0, L0])

    # Solve the system
    solution = solve_ivp(f, [t0, t_final], y0, method='BDF', t_eval=t_eval)
    return solution.y

# Define the objective function
def objective_function(params):
    t_eval = target_viremia_data[:, 0]  # Assume time points in both target datasets are the same
    model_output = model(params, t_eval)

    # Model outputs to compare
    model_viremia = model_output[0]  # Viremia
    model_antibody = model_output[11] + model_output[12]  # Total antibody population

    # Calculate the squared errors
    viremia_error = np.sum((model_viremia - target_viremia_data[:, 1]) ** 2)
    antibody_error = np.sum((model_antibody - target_antibody_data[:, 1]) ** 2)
    
    # Combine the errors
    total_error = viremia_error + antibody_error
    
    return total_error

# Set the bounds for the parameters to be optimized
bounds = [
    (0.1, 10.0),  # Bounds for pi_v
    (1e-5, 1e-1),  # Bounds for alpha_ap
    (1e-3, 1.0),  # Bounds for beta_ap
    (1e-3, 1.0)   # Bounds for delta_apm
]

# Run the differential evolution algorithm
#result = differential_evolution(objective_function, bounds)

# Output the optimal parameters
print('Optimal parameters:', result.x)
print('Objective function value:', result.fun)

# Now you can use these optimized parameters to simulate and save the results
optimal_params = result.x
t_eval = np.linspace(t0, t_final, 150000)
optimized_output = model(optimal_params, t_eval)

# Save the optimized results to a CSV file
filename = "optimized_population_data.csv"
labels = ['Time', 'V', 'Ap', 'Apm', 'Thn', 'The', 'Tkn', 'Tke', 'B', 'Ps', 'Pl', 'Bm', 'IgM', 'IgG', 'L']

with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(labels)  # Write header
    for i in range(len(t_eval)):
        writer.writerow([t_eval[i]] + list(optimized_output[:, i]))  # Write time and populations row by row

print(f"Optimized data successfully saved to {filename}")

