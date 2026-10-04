"""Figura 5: diferencia de los cuadrados de las frecuencias de los modos vs. d^2.

Según la Ec. (10), omega_a^2 - omega_s^2 = (2k/I) d^2, así que la relación es lineal
con pendiente beta = 2k/I. La pendiente ajustada da k_din = I beta / 2.
"""

import matplotlib.pyplot as plt
import numpy as np

# Datos (Tabla III)
d = np.array([0.025, 0.050, 0.076, 0.102])  # m
Ts = np.array([1.26, 1.26, 1.24, 1.24])  # s
Ta = np.array([1.30, 1.19, 1.13, 1.07])  # s

# Incertidumbres instrumentales
dd = np.full_like(d, 0.001)  # m (regla)
dT = np.full_like(Ts, 0.01)  # s (resolución de los periodos)

I, DI = 2.782e-2, 0.014e-2  # kg m^2 (Ec. 7.2)

ws2 = (2 * np.pi / Ts) ** 2
wa2 = (2 * np.pi / Ta) ** 2
y = wa2 - ws2  # s^-2
dy = np.sqrt((2 * ws2 * dT / Ts) ** 2 + (2 * wa2 * dT / Ta) ** 2)

x = d**2  # m^2
dx = 2 * d * dd

# Ajuste lineal por mínimos cuadrados: y = beta x + c
n = len(x)
beta, c = np.polyfit(x, y, 1)
res = y - (beta * x + c)

s = np.sqrt(np.sum(res**2) / (n - 2))
sxx = np.sum((x - x.mean()) ** 2)
dbeta = s / np.sqrt(sxx)
dc = s * np.sqrt(np.sum(x**2) / (n * sxx))
r = np.corrcoef(x, y)[0, 1]

k_din = I * beta / 2
dk_din = 0.5 * np.hypot(beta * DI, I * dbeta)

dy_eff = np.sqrt(dy**2 + (beta * dx) ** 2)

print(f"beta = ({beta:.0f} ± {dbeta:.0f}) s^-2 m^-2")
print(f"c = ({c:.1f} ± {dc:.1f}) s^-2")
print(f"R = {r:.4f}")
print(f"k_din = I beta / 2 = ({k_din:.1f} ± {dk_din:.1f}) N/m")

fig, (ax, axr) = plt.subplots(
    2,
    1,
    figsize=(7, 6.5),
    sharex=True,
    gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08},
)

x_line = np.linspace(0, 0.0115, 200)
ax.plot(
    x_line,
    beta * x_line + c,
    color="C1",
    lw=1.8,
    label=r"Ajuste lineal: $\omega_a^2-\omega_s^2 = \beta\,d^2 + c$",
    zorder=1,
)
ax.errorbar(
    x,
    y,
    xerr=dx,
    yerr=dy,
    fmt="o",
    color="C0",
    ecolor="C0",
    capsize=4,
    ms=6,
    label=r"Datos experimentales ($\pm\delta(d^2)$, $\pm\delta(\omega_a^2-\omega_s^2)$)",
    zorder=2,
)

ax.annotate(
    r"Pendiente $\beta = 2k/I$",
    xy=(0.0075, beta * 0.0075 + c),
    xytext=(0.0040, 8.5),
    fontsize=9.5,
    arrowprops=dict(arrowstyle="->", color="0.3"),
)

ax.set_ylim(-4, 14)
ax.set_ylabel(r"$\omega_a^2-\omega_s^2$ (s$^{-2}$)")
ax.set_title(r"Diferencia de los cuadrados de las frecuencias vs. $d^2$")
ax.grid(alpha=0.3)
ax.legend(loc="upper left", frameon=False, fontsize=9)

axr.axhline(0, color="C1", lw=1.2)
axr.errorbar(
    x,
    res,
    yerr=dy_eff,
    fmt="s",
    color="C2",
    ecolor="C2",
    capsize=4,
    ms=5,
    label="Residuos: dato − ajuste",
)
axr.set_xlabel(r"$d^2$ (m$^2$)")
axr.set_ylabel(r"Residuo (s$^{-2}$)")
axr.grid(alpha=0.3)
axr.legend(loc="upper right", frameon=False, fontsize=8.5)
lim = 1.25 * np.max(np.abs(res) + dy_eff)
axr.set_ylim(-lim, lim)
axr.set_xlim(0, 0.0115)

fig.text(
    0.5,
    0.005,
    "Datos: Tabla III · ω = 2π/T · ajuste por mínimos cuadrados · δT = 0,01 s, δd = 1 mm",
    ha="center",
    fontsize=8,
    color="0.4",
)

fig.savefig("figura5_diferencia_cuadrados.png", dpi=200, bbox_inches="tight")
fig.savefig("figura5_diferencia_cuadrados.pdf", bbox_inches="tight")
