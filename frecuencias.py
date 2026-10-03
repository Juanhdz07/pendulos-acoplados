"""Frecuencias experimentales contra la altura del resorte (Tabla III).

(a) Modos normales: omega_s y omega_a.
(b) Pulsación: omega_p medida, superpuesta con omega_a - omega_s calculada con los
    periodos medidos (deberían coincidir según la Ec. de la pulsación).
"""

import matplotlib.pyplot as plt
import numpy as np

# Datos (Tabla III)
y = np.array([0.365, 0.340, 0.314, 0.288])  # m
Ts = np.array([1.26, 1.26, 1.24, 1.24])  # s
Ta = np.array([1.30, 1.19, 1.13, 1.07])  # s
Tp = np.array([199.0, 49.02, 21.14, 16.1])  # s

# Incertidumbres instrumentales
dy = np.full_like(y, 0.001)  # m (regla)
dT = 0.01  # s (resolución de los periodos)

ws = 2 * np.pi / Ts
wa = 2 * np.pi / Ta
wp = 2 * np.pi / Tp
dws = 2 * np.pi * dT / Ts**2
dwa = 2 * np.pi * dT / Ta**2
dwp = 2 * np.pi * dT / Tp**2

wdiff = wa - ws
dwdiff = np.sqrt(dwa**2 + dws**2)

for i in range(len(y)):
    print(
        f"y = {y[i]:.3f} m: ws = {ws[i]:.3f}, wa = {wa[i]:.3f}, "
        f"wp = {wp[i]:.4f}, wa - ws = {wdiff[i]:.3f} rad/s"
    )

fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4.3))

axa.errorbar(
    y,
    ws,
    xerr=dy,
    yerr=dws,
    fmt="o-",
    color="C0",
    capsize=4,
    ms=6,
    label=r"Modo simétrico $\omega_s$",
)
axa.errorbar(
    y,
    wa,
    xerr=dy,
    yerr=dwa,
    fmt="s-",
    color="C3",
    capsize=4,
    ms=6,
    label=r"Modo antisimétrico $\omega_a$",
)
axa.set_xlabel(r"Altura del resorte $y$ (m)")
axa.set_ylabel(r"Frecuencia angular $\omega$ (rad/s)")
axa.set_title("(a) Modos normales", loc="left")
axa.grid(alpha=0.3)
axa.legend(frameon=False)

axb.errorbar(
    y,
    wp,
    xerr=dy,
    yerr=dwp,
    fmt="o-",
    color="C2",
    capsize=4,
    ms=6,
    label=r"Pulsación medida $\omega_{\mathrm{p}} = 2\pi/T_{\mathrm{p}}$",
)
axb.errorbar(
    y,
    wdiff,
    xerr=dy,
    yerr=dwdiff,
    fmt="^--",
    color="C1",
    capsize=4,
    ms=6,
    label=r"$\omega_a - \omega_s$ (periodos medidos)",
)
axb.axhline(0, color="0.6", lw=0.8)
axb.set_xlabel(r"Altura del resorte $y$ (m)")
axb.set_ylabel(r"Frecuencia angular $\omega$ (rad/s)")
axb.set_title("(b) Pulsación", loc="left")
axb.grid(alpha=0.3)
axb.legend(frameon=False)

fig.text(
    0.5,
    -0.02,
    "Datos: Tabla III · ω = 2π/T · δT = 0,01 s, δy = 1 mm",
    ha="center",
    fontsize=8,
    color="0.4",
)
fig.tight_layout()

fig.savefig("frecuencias.png", dpi=200, bbox_inches="tight")
fig.savefig("frecuencias.pdf", bbox_inches="tight")
