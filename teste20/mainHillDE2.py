# ===========================================
# IMPORTS
# ===========================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, interp1d
from scipy.optimize import differential_evolution

# ===========================================
# 1) PARÂMETROS FIXOS (BASE) — SEU ESTILO
# ===========================================
n       = 70
eq      = 13 + n
t_store = 1000

alpha_l = 2.3
beta_l  = 5.2e-5
delta_l = 1.6e-4

alpha_ap = 1.5e-2
beta_ap  = 9e-1
tau      = 21

c_ap1 = 1
c_ap2 = 40
c_ap3 = 10e6  # será recalculado pelo pico do spline (mantido aqui por clareza)

delta_api = 6.1e-2
gama_api  = 2*3.1e-1
delta_apm = 5.1e-1

alpha_th = 7.17e-2
beta_th  = 1e-6
pi_th    = 1e-7
delta_th = 2.2e0

alpha_tk = 2.17e-1
beta_tk  = 1e-8
pi_tk    = 1e-8
delta_tk = 3e-1

alpha_b = 5.5e1
pi_b1   = 4.83e-7
pi_b2   = 1.27e-6

beta_ps = 4.5e-5
beta_pl = 1e-6
beta_bm = 1e-6

delta_ps = 7.61e1
delta_pl = 3.2e-1
gama_bm  = 4.75e-5

pi_bm1 = 1e-4
pi_bm2 = 2.5e3

c_ps1 = 1.24e2
c_ps2 = 2.3e-1
c_ps3 = 9.03e4

c_pl1 = 5.6e1
c_pl2 = 50.5
c_pl3 = 60.0
c_pl4 = 80.0
c_pl5 = 21.0

delta_IgM = 317.89443611
delta_IgG = 3.39e3

c_apm1 = 1
c_apm2 = 20
c_apm3 = 2e4

# CONDIÇÕES INICIAIS
V0   = 9570.81
L0   = 9.04e6
Ap0  = 0.83e6
Apm0 = 0.0
Thn0 = 1.56e6
The0 = 0.0
Tkn0 = 0.91e6
Tke0 = 0.0
B0   = 1.39e6
Ps0  = 0.0
Pl0  = 0.0
Bm0  = 0.0
A0   = 0.0

# ===========================================
# 2) CARREGAR DADOS EXPERIMENTAIS
# ===========================================
# viremia (para construir spline)
data  = np.genfromtxt('./dados/viremiaPorcoInoculado.csv', delimiter=',', skip_header=1)
told  = data[:, 0]
xold  = data[:, 1]
spline_x = CubicSpline(told, xold)

# alvo IgM para ajuste
target_antibody_data = np.genfromtxt('./dados/plot-IGM.csv', delimiter=',', skip_header=1)
t_exp   = target_antibody_data[:, 0]
IgM_exp = target_antibody_data[:, 1]

# malha temporal da simulação
t0      = 0.0
t_final = 100.0
h       = 0.001
t_eval  = np.linspace(t0, t_final, int((t_final - t0)/h))

# c_ap3 recalculado a partir do pico do spline (mantendo tua lógica)
c_ap3 = float(spline_x(t_eval).max() * 0.8)

# ===========================================
# 3) BASELINE EM DICIONÁRIO (SEM GLOBALS)
# ===========================================
base = dict(
    # dimensões
    eq=eq, t_store=t_store,
    # parâmetros
    alpha_l=alpha_l, beta_l=beta_l, delta_l=delta_l,
    alpha_ap=alpha_ap, beta_ap=beta_ap, c_ap1=c_ap1, c_ap2=c_ap2, c_ap3=c_ap3,
    delta_api=delta_api, gama_api=gama_api, delta_apm=delta_apm,
    alpha_th=alpha_th, beta_th=beta_th, pi_th=pi_th, delta_th=delta_th,
    alpha_tk=alpha_tk, beta_tk=beta_tk, pi_tk=pi_tk, delta_tk=delta_tk,
    alpha_b=alpha_b, pi_b1=pi_b1, pi_b2=pi_b2,
    beta_ps=beta_ps, beta_pl=beta_pl, beta_bm=beta_bm,
    delta_ps=delta_ps, delta_pl=delta_pl, gama_bm=gama_bm,
    pi_bm1=pi_bm1, pi_bm2=pi_bm2,
    c_ps1=c_ps1, c_ps2=c_ps2, c_ps3=c_ps3,
    c_pl1=c_pl1, c_pl2=c_pl2, c_pl3=c_pl3, c_pl4=c_pl4, c_pl5=c_pl5,
    delta_IgM=delta_IgM, delta_IgG=delta_IgG,
    c_apm1=c_apm1, c_apm2=c_apm2, c_apm3=c_apm3,
    # CI
    V0=V0, L0=L0, Ap0=Ap0, Apm0=Apm0, Thn0=Thn0, The0=The0,
    Tkn0=Tkn0, Tke0=Tke0, B0=B0, Ps0=Ps0, Pl0=Pl0, Bm0=Bm0, A0=A0
)

# ===========================================
# 4) PARÂMETROS A AJUSTAR E LIMITES (OS SEUS)
# ===========================================
PARAM_NAMES = ["delta_apm", "c_ps1", "c_ps2", "c_ps3",
               "delta_IgM", "beta_ps", "pi_b2", "alpha_b", "beta_th", "pi_th", "alpha_th"]

delta_apm_bounds = (5e-2, 5e-1)
c_ps1_bounds     = (1.0e1, 6.0e2)
c_ps2_bounds     = (2.3e-2, 2.3e0)
c_ps3_bounds     = (1.03e3, 9.03e4)
delta_IgM_bounds = (3.17e1, 1.17e3)
beta_ps_bounds   = (4.5e-6, 4.5e-4)
pi_b2_bounds     = (1.27e-7, 1.27e-5)
alpha_b_bounds   = (5.5e0, 5.5e2)
beta_th_bounds   = (1e-7, 1e-5)
pi_th_bounds     = (1e-8, 1e-6)
alpha_th_bounds  = (1e-3, 1e-1)


bounds = [
    delta_apm_bounds, c_ps1_bounds, c_ps2_bounds, c_ps3_bounds,
    delta_IgM_bounds, beta_ps_bounds, pi_b2_bounds, alpha_b_bounds, beta_th_bounds, pi_th_bounds, alpha_th_bounds
]

# ===========================================
# 5) CONSTRUTORES: f(t,y) e y0 (dependem de n/tau)
# ===========================================
def make_f(p, d_spline_x):

    def f(t, y):
        dydt = np.zeros(eq)

        # 0) "Vírus" como derivada do spline até dia 54
        dydt[0] = 0.0 if t > 54.0 else float(d_spline_x(t))

        # Saturação para APC (corrigindo c_ap3*c_ap2)
        sat = (p["c_ap1"] * y[0] * p["c_ap2"]) / (p["c_ap3"] * p["c_ap2"] + y[0] * p["c_ap2"])

        # 1) APC naive (Ap)
        dydt[1] = p["alpha_ap"] * (p["Ap0"] - y[1]) - p["beta_ap"] * y[1] * sat

        # 2) "ideal" / pré-cadeia Apm
        dydt[2] = p["beta_ap"] * y[1] * sat - p["delta_apm"] * y[2]

        # 13) Apresentadora Intermediária (início da cadeia)
        dydt[13] = (n * (y[2] - y[13])) / tau

        # Cadeia Apm: y[14]..y[eq-1]
        for i in range(eq - n + 1, eq):
            dydt[i] = (n * (y[i-1] - y[i])) / tau

        Apm_eff = y[eq - 1]  # última célula da cadeia é a efetiva

        # Th naive / efetora
        dydt[3] = p["alpha_th"] * (p["Thn0"] - y[3]) - p["beta_th"] * Apm_eff * y[3]
        dydt[4] = p["beta_th"] * Apm_eff * y[3] + p["pi_th"] * Apm_eff * y[4] - p["delta_th"] * y[4]

        # Tk naive / efetora
        dydt[5] = p["alpha_tk"] * (p["Tkn0"] - y[5]) - p["beta_tk"] * Apm_eff * y[5]
        dydt[6] = p["beta_tk"] * Apm_eff * y[5] + p["pi_tk"] * Apm_eff * y[6] - p["delta_tk"] * y[6]

        # B (com alpha_b e pi_b2 ajustáveis)
        dydt[7] = (
            p["alpha_b"] * (p["B0"] - y[7])
            + p["pi_b2"] * y[4] * y[7]
            - p["beta_ps"] * Apm_eff * y[7]
            - p["beta_pl"] * y[4] * y[7]
            - p["beta_bm"] * y[4] * y[7]
        )

        # Ps, Pl, Bm
        dydt[8]  = p["beta_ps"] * Apm_eff * y[7] - p["delta_ps"] * y[8]
        dydt[9]  = p["beta_pl"] * y[4] * y[7] - p["delta_pl"] * y[9] + p["gama_bm"] * y[10]
        dydt[10] = p["beta_bm"] * y[4] * y[7] + p["pi_bm1"] * y[10] * (1 - y[10] / p["pi_bm2"]) - p["gama_bm"] * y[10]

        # IgM (forma saturante com c_ps1,c_ps2,c_ps3 — equivalente a c_ps1*y8/(c_ps3+y8))
        IgM_sec  = (p["c_ps1"] * y[8] * p["c_ps2"]) / (p["c_ps3"] * p["c_ps2"] + y[8] * p["c_ps2"])
        dydt[11] = IgM_sec - p["delta_IgM"] * y[11]

        # IgG
        dydt[12] = p["c_pl1"] * y[9] - p["delta_IgG"] * y[12]

        return dydt

    return f


def make_y0(p):
    # 13 estados "centrais" + n compartimentos da cadeia Apm
    y0 = np.array([p["V0"], p["Ap0"], p["Apm0"], p["Thn0"], p["The0"], p["Tkn0"], p["Tke0"],
                   p["B0"], p["Ps0"], p["Pl0"], p["Bm0"], p["A0"], p["A0"]], dtype=float)
    for _ in range(n):
        y0 = np.append(y0, p["Apm0"])
    return np.maximum(y0, 1e-12)

# ===========================================
# 6) SIMULADOR E FUNÇÃO DE CUSTO (LOG-SSE IgM)
# ===========================================
def simulate(x):
    # mapear vetor -> dicionário de parâmetros
    p = base.copy()
    mp = dict(zip(PARAM_NAMES, x))

    # tratar n (inteiro >= 1) e tau (>0)

    # copiar parâmetros ajustáveis
    p["delta_apm"] = float(mp["delta_apm"])
    p["c_ps1"]     = float(mp["c_ps1"])
    p["c_ps2"]     = float(mp["c_ps2"])
    p["c_ps3"]     = max(float(mp["c_ps3"]), 1e-9)  # evita zero no denominador
    p["delta_IgM"] = float(mp["delta_IgM"])
    p["beta_ps"]   = float(mp["beta_ps"])
    p["pi_b2"]     = float(mp["pi_b2"])
    p["alpha_b"]   = float(mp["alpha_b"])
    p["beta_th"]   = float(mp["beta_th"])
    p["pi_th"]     = float(mp["pi_th"])
    p["alpha_th"]  = float(mp["alpha_th"])

    # recalcula c_ap3 com base no pico do spline (mantendo sua lógica)
    p["c_ap3"] = float(spline_x(t_eval).max() * 0.8)

    # construir ODE e CI com novo tamanho
    f = make_f(p, spline_x.derivative(1))
    y0_sim = make_y0(p)

    try:
        sol = solve_ivp(
            f, [t0, t_final], y0_sim,
            method='LSODA', t_eval=t_eval,
            rtol=1e-6, atol=1e-9, max_step=1.0
        )
        if (not sol.success) or (sol.y.shape[0] != p["eq"]):
            return None
        Y = np.maximum(sol.y.T, 1e-12)
        if not np.all(np.isfinite(Y)) or np.max(Y) > 1e12:
            return None
        return p, sol.t, Y
    except Exception:
        return None


def cost(x):
    sim = simulate(x)
    if sim is None:
        return 1e12
    p_sim, t_sim, Y = sim
    IgM_sim = Y[:, 11]

    # interpolar no tempo experimental
    interp_IgM = interp1d(t_sim, IgM_sim, kind='linear', fill_value='extrapolate', bounds_error=False)
    IgM_model  = np.maximum(interp_IgM(t_exp), 1e-12)
    IgM_obs    = np.maximum(IgM_exp, 1e-12)

    # erro em log (coerente com escala)
    err = np.log(IgM_obs) - np.log(IgM_model)
    SSE = float(np.sum(err**2))

    # penalizações adicionais (solução absurda)
    if np.any(IgM_model > 1e9):
        SSE += 1e10

    return SSE

# ===========================================
# 7) AJUSTE COM DE
# ===========================================
result = differential_evolution(
    cost, bounds,
    strategy='best2bin',
    seed=42,
    tol=1e-7,
    mutation=(0.6, 1.6),
    recombination=0.9,
    maxiter=2000, 
    popsize=100,
    polish=True,
    disp=True,
    workers=-1,
)

print("\n===== RESULTADO DO AJUSTE =====")
print(f"SSE (log-IgM) = {result.fun:.6g}")
for nm, val in zip(PARAM_NAMES, result.x):
    print(f"{nm:>10s} = {val:.6g}")

# ===========================================
# 8) SIMULA MELHOR SOLUÇÃO, SALVA CSV E PLOTA
# ===========================================
best = result.x
sim_best = simulate(best)
if sim_best is None:
    raise RuntimeError("Falha ao simular com os melhores parâmetros.")

p_best, t_values, Y_values = sim_best

# salvar CSV do melhor (com rótulos consistentes com n ótimo)
filename = "output_python.csv"
labels = ['Time', 'V', 'Ap','Api', 'Apm', 'Thn', 'The', 'Tkn', 'Tke', 'B', 'Ps', 'Pl', 'Bm', 'IgM', 'IgG']
for i in range(n):
    labels.append('Ap'+str(i))

with open(filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(labels)
    for i in range(len(t_values)):
        row = [t_values[i]] + list(Y_values[i, :14])  # 0..13 fixos
        # anexar cadeia Apm (n compartimentos)
        row += list(Y_values[i, 14:14 + n])
        writer.writerow(row)

# plot IgM exp vs sim (log, coerente com custo)
plt.figure()
plt.plot(t_exp, IgM_exp, 'o', label='IgM exp')
plt.plot(t_values, Y_values[:, 11], '-', label='IgM sim (best)')
plt.yscale('log')
plt.xlabel('Dias')
plt.ylabel('IgM')
plt.legend(loc='best')
plt.title('Ajuste IgM (escala log)')
plt.tight_layout()
plt.show()

print(f"\nDados salvos em: {filename}")

