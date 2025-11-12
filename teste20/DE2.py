import numpy as np
from scipy.optimize import differential_evolution
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d

# ============================================================
# Base model definitions
# ============================================================

# Initial constants
n = 70
eq = 13 + n
t_store = 1000

alpha_l = 2.3
beta_l = 5.2e-5
delta_l = 1.6e-4
alpha_ap = 1.5e-2
beta_ap = 9e-1
tau = 21
c_ap1 = 1
c_ap2 = 40
c_ap3 = 10e6
delta_api = 6.1e-2
gama_api = 2 * 3.1e-1
delta_apm = 5.1e-1
alpha_th = 7.17e-2
beta_th = 1e-6
pi_th = 1e-7
delta_th = 1e-1
alpha_tk = 5.17e-2
beta_tk = 1e-7
pi_tk = 1e-7
delta_tk = 1e-1
alpha_b = 4.1e1
pi_b2 = 1.27e-6
beta_ps = 4.5e-5
beta_pl = 4.5e-5
beta_bm = 4.5e-5
delta_ps = 3.4e1
delta_pl = 3.4e1
delta_bm = 3.4e1
gama_bm = 4.1e-2
pi_bm1 = 2.3e-4
pi_bm2 = 1.7e1
delta_IgM = 1e2
c_ps1 = 5e2
c_ps2 = 1
c_ps3 = 1
c_pl1 = 3e2
delta_IgG = 5e1

# Initial conditions
V0 = 1e6
Ap0 = 1e6
Apm0 = 1e6
Thn0 = 1e6
The0 = 1e6
Tkn0 = 1e6
Tke0 = 1e6
B0 = 1e6
Ps0 = 1e6
Pl0 = 1e6
Bm0 = 1e6
A0 = 0

# Load viremia spline
data = np.genfromtxt('./dados/viremiaPorcoInoculado.csv', delimiter=',', skip_header=1)
told = data[:, 0]
xold = data[:, 1]
spline_x = CubicSpline(told, xold)

# Experimental IgM data
IgM_exp_data = np.genfromtxt('./dados/plot-IGM.csv', delimiter=',', skip_header=1)
t_exp = IgM_exp_data[:, 0]
IgM_exp = IgM_exp_data[:, 1]

# Time range
t0 = 0
t_final = 100
t_eval = np.linspace(t0, t_final, 1000)

# ============================================================
# Model equations (parameterized)
# ============================================================

def f(t, y, d_spline_x, params=None):
    if params is None:
        raise ValueError("params dictionary must be provided")

    n = int(params["n"])
    eq = params["eq"]
    #tau = params["tau"]
    delta_apm = params["delta_apm"]
    c_ps1 = params["c_ps1"]
    #c_ps2 = params["c_ps2"]
    #c_ps3 = params["c_ps3"]
    delta_IgM = params["delta_IgM"]
    beta_ps = params["beta_ps"]
    pi_b2 = params["pi_b2"]
    alpha_b = params["alpha_b"]
    beta_th = params["beta_th"]
    pi_th = params["pi_th"]
    alpha_th = params["alpha_th"]
    #c_ap3 = params["c_ap3"]

    dydt = np.zeros(eq)

    # Virus
    dydt[0] = 0 if t > 54 else d_spline_x(t)
    # APC naive
    dydt[1] = alpha_ap * (Ap0 - y[1]) - beta_ap * y[1] * ((c_ap1 * y[0] ** c_ap2) / (c_ap3 ** c_ap2 + y[0] ** c_ap2))
    # APC ideal
    dydt[2] = beta_ap * y[1] * ((c_ap1 * y[0] ** c_ap2) / (c_ap3 ** c_ap2 + y[0] ** c_ap2)) - delta_apm * y[2]
    # APC intermediate
    dydt[13] = (n * (y[2] - y[13])) / tau
    # APC mature
    for i in range(eq - n + 1, eq):
        dydt[i] = (n * (y[i - 1] - y[i])) / tau
    # Th naive
    dydt[3] = alpha_th * (Thn0 - y[3]) - beta_th * y[eq - 1] * y[3]
    # Th effector
    dydt[4] = beta_th * y[eq - 1] * y[3] + pi_th * y[eq - 1] * y[4] - delta_th * y[4]
    # B
    dydt[7] = alpha_b * (B0 - y[7]) + pi_b2 * y[4] * y[7] - beta_ps * y[eq - 1] * y[7]
    # Ps
    dydt[8] = beta_ps * y[eq - 1] * y[7] - delta_ps * y[8]
    # IgM
    dydt[11] = (c_ps1 * (y[8] ** c_ps2) / (c_ps3 ** c_ps2 + y[8] ** c_ps2)) - delta_IgM * y[11]

    return dydt

# ============================================================
# Error metric
# ============================================================

def erro(t_values, y_values):
    IgM_model = interp1d(t_values, y_values[:, 11], kind="linear", fill_value="extrapolate")(t_exp)
    l2_rel_error = np.linalg.norm(IgM_model - IgM_exp) / np.linalg.norm(IgM_exp)
    return l2_rel_error * 100.0

# ============================================================
# Objective function (no globals)
# ============================================================

def funcaoObjetivo(p):
    keys = ["n", "tau", "delta_apm", "c_ps1", "c_ps2", "c_ps3",
            "delta_IgM", "beta_ps", "pi_b2", "alpha_b",
            "beta_th", "pi_th", "alpha_th"]
    params = {k: float(p[i]) for i, k in enumerate(keys)}
    params["n"] = max(1, int(round(params["n"])))
    params["eq"] = 13 + params["n"]
    params["c_ap3"] = float(spline_x(t_eval).max() * 0.8)

    y0_local = np.array([V0, Ap0, Apm0, Thn0, The0, Tkn0, Tke0,
                         B0, Ps0, Pl0, Bm0, A0, A0])
    for _ in range(13, params["eq"]):
        y0_local = np.append(y0_local, Apm0)

    try:
        sol = solve_ivp(f, [t0, t_final], y0_local, method="LSODA", t_eval=t_eval,
                        args=(spline_x.derivative(1), params),
                        rtol=1e-6, atol=1e-9)
        if not sol.success:
            return 1e12
    except Exception:
        return 1e12

    try:
        error = erro(sol.t, sol.y.T)
    except Exception:
        return 1e12

    return error if np.isfinite(error) else 1e12

# ============================================================
# Parameter bounds & DE optimization
# ============================================================

PARAM_NAMES = ["n", "tau", "delta_apm", "c_ps1", "c_ps2", "c_ps3",
               "delta_IgM", "beta_ps", "pi_b2", "alpha_b",
               "beta_th", "pi_th", "alpha_th"]

bounds = [
    (0, 50), (0, 50), (5e-2, 5e0), (1.0e1, 1.0e3), (0, 50),
    (0, 240000), (3.17e1, 1.17e3), (4.5e-6, 4.5e-4), (1.27e-7, 1.27e-5),
    (5.5e0, 5.5e2), (1e-7, 1e-5), (1e-8, 1e-6), (1e-3, 1e-1)
]

result = differential_evolution(
    funcaoObjetivo,
    bounds,
    strategy="best1bin",
    popsize=15,
    mutation=(0.5, 1.0),
    recombination=0.7,
    seed=42,
    tol=1e-6,
    maxiter=200,
    polish=True,
    updating="deferred",
    workers=-1
)

print("\n=== RESULTADOS DE ===")
print("Sucesso:", result.success)
print("Mensagem:", result.message)
print("Erro (objetivo) mínimo:", result.fun)
best_p = result.x
for name, val in zip(PARAM_NAMES, best_p):
    print(f"{name}* = {val:.6g}")

# ============================================================
# Re-simulate best parameters
# ============================================================

params_best = {k: float(best_p[i]) for i, k in enumerate(PARAM_NAMES)}
params_best["n"] = max(1, int(round(params_best["n"])))
params_best["eq"] = 13 + params_best["n"]
params_best["c_ap3"] = float(spline_x(t_eval).max() * 0.8)

y0_best = np.array([V0, Ap0, Apm0, Thn0, The0, Tkn0, Tke0,
                    B0, Ps0, Pl0, Bm0, A0, A0])
for _ in range(13, params_best["eq"]):
    y0_best = np.append(y0_best, Apm0)

sol_best = solve_ivp(f, [t0, t_final], y0_best, method="LSODA",
                     t_eval=t_eval, args=(spline_x.derivative(1), params_best),
                     rtol=1e-6, atol=1e-9)

IgM_best = np.interp(t_exp, sol_best.t, sol_best.y.T[:, 11])

plt.figure()
plt.plot(sol_best.t, sol_best.y.T[:, 11], label="IgM (modelo)")
plt.scatter(t_exp, IgM_exp, s=30, label="IgM (dados)")
plt.xlabel("Tempo (dias)")
plt.ylabel("IgM")
plt.legend()
plt.title("Ajuste de parâmetros via DE (IgM)")
plt.show()

