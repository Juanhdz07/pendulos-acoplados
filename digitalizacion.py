"""Digitalización de la curva ángulo vs. tiempo de la captura de Logger Pro.

Para cada columna de píxeles de la curva se toma el centro del trazo como el ángulo
y la mitad de su grosor como incertidumbre de lectura.
"""

import numpy as np
from PIL import Image

# Calibración de ejes en la captura (píxel -> valor)
T_PX = np.array([149, 470, 791])
T_VAL = np.array([0.0, 10.0, 20.0])
Y_PX = np.array([404, 466, 526, 587])
Y_VAL = np.array([80.0, 60.0, 40.0, 20.0])

# Región de la curva dentro del área de la gráfica
ROW_MIN, ROW_MAX = 406, 622
COL_MIN, COL_MAX = 152, 560
DARK = 110  # nivel de gris máximo del trazo (fondo y cuadrícula son más claros)


def digitalizar(path="logger_pro_pulsacion.png", step=1):
    """Devuelve (t, theta, dtheta) en s y grados, tomando una columna de cada `step`."""
    t_slope, t_off = np.polyfit(T_PX, T_VAL, 1)
    y_slope, y_off = np.polyfit(Y_PX, Y_VAL, 1)
    img = np.asarray(Image.open(path).convert("L"), dtype=float)

    t, th, dth = [], [], []
    for col in range(COL_MIN, COL_MAX, step):
        if np.min(np.abs(T_PX - col)) <= 2:  # líneas verticales de la cuadrícula
            continue
        rows = np.nonzero(img[ROW_MIN:ROW_MAX, col] < DARK)[0] + ROW_MIN
        if rows.size == 0:
            continue
        center = 0.5 * (rows.min() + rows.max())
        half = max(0.5 * (rows.max() - rows.min()), 0.5)
        t.append(t_slope * col + t_off)
        th.append(y_slope * center + y_off)
        dth.append(abs(y_slope) * half)

    return np.array(t), np.array(th), np.array(dth)
