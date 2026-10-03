"""Ajuste lineal m(x) para la constante del resorte, con barras de error y residuos.

Tabla II: alargamiento x del resorte en función de la masa colgada m.
En equilibrio k x = m g, por lo que la pendiente b de m(x) da k = b g.
"""

import matplotlib.pyplot as plt
import numpy as np

# Datos (Tabla II)
x = np.array([0.000, 0.026, 0.029, 0.032, 0.054])  # m
m = np.array([0.00, 60.00, 70.00, 80.00, 130.00])  # g

# Incertidumbres instrumentales: regla (1 mm) y balanza (0,01 g)
dx = np.full_like(x, 0.001)  # m
dm = np.full_like(m, 0.01)  # g

G = 9.8  # m/s^2

# Ajuste lineal por mínimos cuadrados: m = b x + c
n = len(x)
b, c = np.polyfit(x, m, 1)
m_fit = b * x + c
res = m - m_fit

s = np.sqrt(np.sum(res**2) / (n - 2))
sxx = np.sum((x - x.mean()) ** 2)
db = s / np.sqrt(sxx)
dc = s * np.sqrt(np.sum(x**2) / (n * sxx))
r = np.corrcoef(x, m)[0, 1]

k = b / 1000 * G  # N/m
dk = db / 1000 * G

# Incertidumbre efectiva en m: incluye la de x proyectada con la pendiente
dm_eff = np.sqrt(dm**2 + (b * dx) ** 2)

print(f"b = ({b:.0f} ± {db:.0f}) g/m")
print(f"c = ({c:.1f} ± {dc:.1f}) g")
print(f"R = {r:.4f}")
print(f"k = ({k:.1f} ± {dk:.1f}) N/m")

fig, (ax, axr) = plt.subplots(
    2,
    1,
    figsize=(7, 6.5),
    sharex=True,
    gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08},
)

x_line = np.linspace(-0.002, 0.058, 200)
ax.plot(
    x_line,
    b * x_line + c,
    color="C1",
    lw=1.8,
    label=rf"Ajuste lineal: $m = b\,x + c$",
    zorder=1,
)
ax.errorbar(
    x,
    m,
    xerr=dx,
    yerr=dm,
    fmt="o",
    color="C0",
    ecolor="C0",
    capsize=4,
    ms=6,
    label=r"Datos experimentales ($\pm\delta x$, $\pm\delta m$)",
    zorder=2,
)

ax.set_ylabel(r"Masa colgada $m$ (g)")
ax.set_title("Constante del resorte: masa colgada vs. alargamiento")
ax.grid(alpha=0.3)
ax.legend(loc="upper left", frameon=False)

ax.text(
    0.97,
    0.05,
    (
        rf"$b = ({b:.0f} \pm {db:.0f})$ g/m" + "\n"
        rf"$c = ({c:.1f} \pm {dc:.1f})$ g" + "\n"
        rf"$R = {r:.4f}$" + "\n"
        rf"$k = b\,g = ({k:.1f} \pm {dk:.1f})$ N/m"
    ),
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=10,
    bbox=dict(boxstyle="round", fc="white", ec="0.7"),
)

ax.annotate(
    "Barras de error\n(δx = 1 mm)",
    xy=(x[1] - dx[1], m[1]),
    xytext=(0.004, 72),
    fontsize=9,
    ha="left",
    arrowprops=dict(arrowstyle="->", color="0.3"),
)
ax.annotate(
    "Recta de ajuste",
    xy=(0.045, b * 0.045 + c),
    xytext=(0.047, 80),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="0.3"),
)

axr.axhline(0, color="C1", lw=1.2)
axr.errorbar(
    x,
    res,
    yerr=dm_eff,
    fmt="s",
    color="C2",
    ecolor="C2",
    capsize=4,
    ms=5,
    label=r"Residuos $m - m_{\mathrm{ajuste}}$ ($\pm\sqrt{\delta m^2 + (b\,\delta x)^2}$)",
)
axr.set_xlabel(r"Alargamiento $x$ (m)")
axr.set_ylabel("Residuo (g)")
axr.grid(alpha=0.3)
axr.legend(loc="upper left", frameon=False, fontsize=8.5)
lim = 1.25 * np.max(np.abs(res) + dm_eff)
axr.set_ylim(-lim, lim)
axr.set_xlim(-0.003, 0.058)

fig.text(
    0.5,
    0.005,
    "Datos: Tabla II · ajuste por mínimos cuadrados · g = 9,8 m/s²",
    ha="center",
    fontsize=8,
    color="0.4",
)

fig.savefig("constante_resorte.png", dpi=200, bbox_inches="tight")
fig.savefig("constante_resorte.pdf", bbox_inches="tight")
