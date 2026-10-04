"""Ángulo vs. tiempo de un péndulo con pulsación (resorte en y = 0,288 m).

Los datos se digitalizan de la captura de Logger Pro (logger_pro_pulsacion.png):
para cada columna de píxeles de la curva se toma el centro del trazo como el ángulo
y la mitad de su grosor como incertidumbre de lectura.

Se ajusta el modelo de pulsación
    theta(t) = theta0 + A cos(Omega t + phi1) cos(omega t + phi2),
con omega = (omega_a + omega_s)/2 y Omega = (omega_a - omega_s)/2.
"""

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.optimize import curve_fit

# Calibración de ejes en la captura (píxel -> valor)
T_PX = np.array([149, 470, 791])
T_VAL = np.array([0.0, 10.0, 20.0])
Y_PX = np.array([404, 466, 526, 587])
Y_VAL = np.array([80.0, 60.0, 40.0, 20.0])
t_slope, t_off = np.polyfit(T_PX, T_VAL, 1)
y_slope, y_off = np.polyfit(Y_PX, Y_VAL, 1)

# Región de la curva dentro del área de la gráfica
ROW_MIN, ROW_MAX = 406, 622
COL_MIN, COL_MAX = 152, 560
DARK = 110  # nivel de gris máximo del trazo (fondo y cuadrícula son más claros)
STEP = 3  # tomar una columna de cada STEP

img = np.asarray(Image.open("logger_pro_pulsacion.png").convert("L"), dtype=float)

t_list, th_list, dth_list = [], [], []
for col in range(COL_MIN, COL_MAX, STEP):
    if np.min(np.abs(T_PX - col)) <= 2:  # líneas verticales de la cuadrícula
        continue
    rows = np.nonzero(img[ROW_MIN:ROW_MAX, col] < DARK)[0] + ROW_MIN
    if rows.size == 0:
        continue
    center = 0.5 * (rows.min() + rows.max())
    half = max(0.5 * (rows.max() - rows.min()), 0.5)
    t_list.append(t_slope * col + t_off)
    th_list.append(y_slope * center + y_off)
    dth_list.append(abs(y_slope) * half)

t = np.array(t_list)
th = np.array(th_list)
dth = np.array(dth_list)


def modelo(t, theta0, A, Omega, phi1, omega, phi2):
    return theta0 + A * np.cos(Omega * t + phi1) * np.cos(omega * t + phi2)


# Valores iniciales: envolvente con mínimos en t ~ 1 s y t ~ 9 s, periodo rápido ~ 1,15 s
Omega0 = np.pi / 8.05
p0 = [th.mean(), 8.0, Omega0, np.pi / 2 - Omega0 * 1.0, 2 * np.pi / 1.15, 0.0]
popt, pcov = curve_fit(modelo, t, th, p0=p0, sigma=dth, absolute_sigma=True, maxfev=20000)
perr = np.sqrt(np.diag(pcov))
theta0, A, Omega, phi1, omega, phi2 = popt
dOmega, domega = perr[2], perr[4]

res = th - modelo(t, *popt)
chi2_red = np.sum((res / dth) ** 2) / (len(t) - len(popt))

ws = omega - abs(Omega)
wa = omega + abs(Omega)
dw = np.hypot(domega, dOmega)
T_env = np.pi / abs(Omega)
dT_env = T_env * dOmega / abs(Omega)

print(f"Puntos digitalizados: {len(t)}")
print(f"theta0 = {theta0:.2f} ± {perr[0]:.2f} grados, A = {abs(A):.2f} ± {perr[1]:.2f} grados")
print(f"omega (portadora) = {omega:.3f} ± {domega:.3f} rad/s")
print(f"Omega (envolvente) = {abs(Omega):.4f} ± {dOmega:.4f} rad/s")
print(f"omega_s = {ws:.3f} ± {dw:.3f} rad/s, omega_a = {wa:.3f} ± {dw:.3f} rad/s")
print(f"Periodo de la envolvente pi/Omega = {T_env:.2f} ± {dT_env:.2f} s")
print(f"chi2 reducido = {chi2_red:.2f}")

fig, (ax, axr) = plt.subplots(
    2,
    1,
    figsize=(8, 6.5),
    sharex=True,
    gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08},
)

t_line = np.linspace(t.min(), t.max(), 1500)
ax.plot(
    t_line,
    modelo(t_line, *popt),
    color="C1",
    lw=1.6,
    label=r"Regresión: $\theta_0 + A\cos(\Omega t+\varphi_1)\cos(\omega t+\varphi_2)$",
    zorder=1,
)
env = abs(A) * np.abs(np.cos(Omega * t_line + phi1))
ax.plot(t_line, theta0 + env, color="0.5", lw=1, ls="--", label="Envolvente", zorder=1)
ax.plot(t_line, theta0 - env, color="0.5", lw=1, ls="--", zorder=1)
ax.errorbar(
    t,
    th,
    yerr=dth,
    fmt="o",
    color="C0",
    ecolor="C0",
    elinewidth=0.8,
    capsize=2,
    ms=3,
    label=r"Datos (digitalizados de Logger Pro, $\pm\delta\theta$)",
    zorder=2,
)

ax.annotate(
    "Intercambio de energía:\nla amplitud crece y decae",
    xy=(4.3, theta0 + abs(A) * 0.98),
    xytext=(6.6, 77),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="0.3"),
)

ax.set_ylabel(r"Ángulo $\theta$ (°)")
ax.set_title(r"Pulsación: ángulo de un péndulo vs. tiempo ($y = 0{,}288$ m)")
ax.set_ylim(38, 82)
ax.grid(alpha=0.3)
ax.legend(loc="lower left", frameon=False, fontsize=8.5)

axr.axhline(0, color="C1", lw=1.2)
axr.errorbar(
    t,
    res,
    yerr=dth,
    fmt="s",
    color="C2",
    ecolor="C2",
    elinewidth=0.8,
    capsize=2,
    ms=3,
    label="Residuos: dato − regresión",
)
axr.set_xlabel(r"Tiempo $t$ (s)")
axr.set_ylabel("Residuo (°)")
axr.grid(alpha=0.3)
axr.legend(loc="upper right", frameon=False, fontsize=8.5)
lim = 1.2 * np.max(np.abs(res) + dth)
axr.set_ylim(-lim, lim)

fig.text(
    0.5,
    0.005,
    "Datos digitalizados de la gráfica de Logger Pro · δθ = medio grosor del trazo · ajuste por mínimos cuadrados ponderados",
    ha="center",
    fontsize=8,
    color="0.4",
)

fig.savefig("pulsacion.png", dpi=200, bbox_inches="tight")
fig.savefig("pulsacion.pdf", bbox_inches="tight")
