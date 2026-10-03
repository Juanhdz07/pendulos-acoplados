"""Modos normales de los péndulos acoplados (Tabla III), con barras de error y residuos.

Para dos péndulos físicos idénticos acoplados por un resorte a distancia d del eje,
omega_a^2 = omega_s^2 + 2 k d^2 / I, así que omega_a^2 - omega_s^2 es lineal en d^2.
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
ws2 = (2 * np.pi / Ts) ** 2
wa2 = (2 * np.pi / Ta) ** 2
y = wa2 - ws2  # s^-2
dy = np.sqrt((2 * ws2 * dT / Ts) ** 2 + (2 * wa2 * dT / Ta) ** 2)

x = d**2  # m^2
dx = 2 * d * dd

# Ajuste lineal por mínimos cuadrados: y = b x + c
n = len(x)
b, c = np.polyfit(x, y, 1)
y_fit = b * x + c
res = y - y_fit

s = np.sqrt(np.sum(res**2) / (n - 2))
sxx = np.sum((x - x.mean()) ** 2)
db = s / np.sqrt(sxx)
dc = s * np.sqrt(np.sum(x**2) / (n * sxx))
r = np.corrcoef(x, y)[0, 1]

# Incertidumbre efectiva en y: incluye la de x proyectada con la pendiente
dy_eff = np.sqrt(dy**2 + (b * dx) ** 2)

print(f"b = ({b:.0f} ± {db:.0f}) s^-2 m^-2")
print(f"c = ({c:.1f} ± {dc:.1f}) s^-2")
print(f"R = {r:.4f}")

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
    b * x_line + c,
    color="C1",
    lw=1.8,
    label=r"Ajuste lineal: $\omega_a^2-\omega_s^2 = b\,d^2 + c$",
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
for i in range(n):
    ax.annotate(
        f"{i + 1}",
        (x[i], y[i]),
        xytext=(6, -12),
        textcoords="offset points",
        fontsize=9,
        color="C0",
    )

ax.set_ylim(-4, 14)
ax.set_ylabel(r"$\omega_a^2-\omega_s^2$ (s$^{-2}$)")
ax.set_title("Modos normales: separación de frecuencias vs. distancia del resorte")
ax.grid(alpha=0.3)
ax.legend(loc="upper left", frameon=False, fontsize=9)

ax.annotate(
    "Punto 1: $T_a > T_s$\n(acoplamiento débil)",
    xy=(x[0], y[0]),
    xytext=(0.0016, -3.3),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="0.3"),
)
ax.annotate(
    "Recta de ajuste",
    xy=(0.0075, b * 0.0075 + c),
    xytext=(0.0045, 8.5),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="0.3"),
)

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

fig.savefig("modos_normales.png", dpi=200, bbox_inches="tight")
fig.savefig("modos_normales.pdf", bbox_inches="tight")
