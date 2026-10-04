"""Transformada de Fourier de la señal de pulsación (resorte en y = 0,288 m).

La señal digitalizada de Logger Pro se interpola a un paso de tiempo uniforme, se le
resta el valor medio y se calcula el espectro de amplitud con la FFT. Una pulsación
es la suma de los dos modos normales, así que el espectro debe mostrar dos picos:
uno en omega_s y otro en omega_a.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from digitalizacion import digitalizar

t, th, _ = digitalizar(step=1)

dt = np.median(np.diff(t))
t_u = np.arange(t.min(), t.max(), dt)
th_u = np.interp(t_u, t, th)
x = th_u - th_u.mean()

n = len(x)
T_total = n * dt
N_PAD = 16 * n  # relleno con ceros: interpola el espectro, no mejora la resolución

spec = np.abs(np.fft.rfft(x, N_PAD)) * 2 / n  # amplitud (grados)
f = np.fft.rfftfreq(N_PAD, dt)  # Hz
w = 2 * np.pi * f  # rad/s

banda = (w > 3) & (w < 8)
picos, _ = find_peaks(spec * banda, height=0.3 * spec[banda].max())
picos = picos[np.argsort(spec[picos])[::-1][:2]]
picos = np.sort(picos)

dw = 2 * np.pi / T_total  # resolución en frecuencia angular

print(f"Duración de la señal: {T_total:.2f} s, dt = {dt * 1000:.1f} ms, {n} puntos")
print(f"Resolución: Δω = 2π/T = {dw:.3f} rad/s")
for p in picos:
    print(f"Pico: ω = {w[p]:.3f} rad/s (f = {f[p]:.3f} Hz), amplitud = {spec[p]:.2f} grados")
if len(picos) == 2:
    print(f"Separación: Δω = {w[picos[1]] - w[picos[0]]:.3f} rad/s")

fig, ax = plt.subplots(figsize=(7.5, 4.3))
ax.plot(w, spec, color="C0", lw=1.6, label="Espectro de amplitud |FFT|")

nombres = [r"$\omega_s$", r"$\omega_a$"]
for nombre, p in zip(nombres, picos):
    ax.plot(w[p], spec[p], "o", color="C3", ms=6)
    ax.annotate(
        f"{nombre} = {w[p]:.2f} rad/s",
        xy=(w[p], spec[p]),
        xytext=(0, 10),
        textcoords="offset points",
        ha="center",
        fontsize=9.5,
    )

ax.set_xlim(0, 12)
ax.set_ylim(0, 1.25 * spec[banda].max())
ax.set_xlabel(r"Frecuencia angular $\omega$ (rad/s)")
ax.set_ylabel("Amplitud (°)")
ax.set_title(
    "Transformada de Fourier de la pulsación\n"
    r"Pulso de la columna T(puls) de la tabla: $T_{\mathrm{puls}} = 16{,}1$ s ($y = 0{,}288$ m)",
    fontsize=11,
)
ax.grid(alpha=0.3)
ax.legend(loc="upper right", frameon=False)

fig.text(
    0.5,
    -0.02,
    f"Señal digitalizada de Logger Pro · {T_total:.1f} s · resolución Δω = 2π/T = {dw:.2f} rad/s · media restada",
    ha="center",
    fontsize=8,
    color="0.4",
)
fig.tight_layout()

fig.savefig("fourier.png", dpi=200, bbox_inches="tight")
fig.savefig("fourier.pdf", bbox_inches="tight")
