import numpy as np

mu = 398600     # km^3/s^2

# Conversión de elementos orbitales (COE) a vectores de posición y velocidad en ECI.
def COE2RV(p: float, e: float, i: float, raan: float, w: float, nu: float,
           u: float = None, lambda_true: float = None, w_true: float = None):

    if e == 0 and i == 0:
        w = 0
        raan = 0
        nu = lambda_true
    if e == 0 and i != 0:
        w = 0
        nu = u
    if e != 0 and i == 0:
        raan = 0
        w = w_true

    r_pqw = [
        p * np.cos(nu) / (1 + e * np.cos(nu)),
        p * np.sin(nu) / (1 + e * np.cos(nu)),
        0
    ]

    v_pqw = [
        -np.sqrt(mu / p) * np.sin(nu),
        np.sqrt(mu / p) * (e + np.cos(nu)),
        0
    ]

    ijk_pqw = np.array([
        [
            np.cos(raan) * np.cos(w) - np.sin(raan) * np.sin(w) * np.cos(i),
            -np.cos(raan) * np.sin(w) - np.sin(raan) * np.cos(w) * np.cos(i),
            np.sin(raan) * np.sin(i)
        ],
        [
            np.sin(raan) * np.cos(w) + np.cos(raan) * np.sin(w) * np.cos(i),
            -np.sin(raan) * np.sin(w) + np.cos(raan) * np.cos(w) * np.cos(i),
            -np.cos(raan) * np.sin(i)
        ],
        [
            np.sin(w) * np.sin(i),
            np.cos(w) * np.sin(i),
            np.cos(i)
        ]
    ])

    r_ijk = np.dot(ijk_pqw, r_pqw)
    v_ijk = np.dot(ijk_pqw, v_pqw)

    return np.round(r_ijk, 5), np.round(v_ijk, 5)

# Datos de Prueba

# Elementos orbitales de ejemplo
p = 11067.790           # km
e = 0.83285
i = 87.87               # grados
raan = 227.89           # grados
w = 53.38               # grados
nu = 92.335             # grados

# Convertir ángulos a radianes
i = np.deg2rad(i)
raan = np.deg2rad(raan)
w = np.deg2rad(w)
nu = np.deg2rad(nu)

# Ejecutar conversión
r_vec, v_vec = COE2RV(p, e, i, raan, w, nu)

print("> Vector de posición r_ijk (km): ", r_vec)
print("> Vector de velocidad v_ijk (km/s): ", v_vec)
