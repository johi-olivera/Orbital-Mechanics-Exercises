import numpy as np

mu = 398600     # km^3/s^2

# Calcula los Elementos Orbitales Clásicos (COE) a partir del vector de posición y velocidad.
def convertir_rv_a_coe(r_ijk, v_ijk):
    """
    Parámetros:
    - r_ijk: vector de posición en ECI [km]
    - v_ijk: vector de velocidad en ECI [km/s]

    Retorna:
    - p: semilato recto [km]
    - a: semieje mayor [km]
    - e: excentricidad
    - i: inclinación [rad]
    - raan: ascensión recta del nodo ascendente [rad]
    - w: argumento del perigeo [rad]
    - nu: anomalía verdadera [rad]
    - w_true: argumento del perigeo verdadero [rad] (solo si órbita elíptica ecuatorial)
    - u: argumento de la latitud [rad] (solo si órbita circular inclinada)
    - lambda_true: longitud verdadera [rad] (solo si órbita circular ecuatorial)
    """

    # Conversión a arrays NumPy
    r_ijk = np.array(r_ijk, dtype=float)
    v_ijk = np.array(v_ijk, dtype=float)

    # Normas de vectores
    r = np.linalg.norm(r_ijk)
    v = np.linalg.norm(v_ijk)

    # Vector momento angular específico
    h_vec = np.cross(r_ijk, v_ijk)
    h = np.linalg.norm(h_vec)

    # Vector nodo ascendente
    n_vec = np.cross([0, 0, 1], h_vec)
    n = np.linalg.norm(n_vec)

    # Vector de excentricidad
    e_vec = ((v**2 - mu / r) * r_ijk - np.dot(r_ijk, v_ijk) * v_ijk) / mu
    e = np.linalg.norm(e_vec)

    # Energía específica orbital
    energia = v**2 / 2 - mu / r

    # Semieje mayor y semilato recto
    if e != 1.0:
        a = -mu / (2 * energia)
        p = a * (1 - e**2)
    else:
        a = np.inf
        p = h**2 / mu

    # Inclinación orbital
    i = np.arccos(h_vec[2] / h)

    # Ascensión recta del nodo ascendente (RAAN)
    raan = np.arccos(n_vec[0] / n)
    if n_vec[1] < 0:
        raan = 2 * np.pi - raan

    # Argumento del perigeo
    w = np.arccos(np.dot(n_vec, e_vec) / (n * e))
    if e_vec[2] < 0:
        w = 2 * np.pi - w

    # Anomalía verdadera
    nu = np.arccos(np.dot(e_vec, r_ijk) / (e * r))
    if np.dot(r_ijk, v_ijk) < 0:
        nu = 2 * np.pi - nu

    # Inicialización de casos especiales
    w_true = None
    u = None
    lambda_true = None

    # Casos especiales
    if e < 1 and i == 0:
        # Órbita elíptica ecuatorial
        w_true = np.arccos(e_vec[0] / e)
        if e_vec[1] < 0:
            w_true = 2 * np.pi - w_true

    elif e == 0 and i != 0:
        # Órbita circular inclinada
        u = np.arccos(np.dot(n_vec, r_ijk) / (n * r))
        if r_ijk[2] < 0:
            u = 2 * np.pi - u

    elif e == 0 and i == 0:
        # Órbita circular ecuatorial
        lambda_true = np.arccos(r_ijk[0] / r)
        if r_ijk[1] < 0:
            lambda_true = 2 * np.pi - lambda_true

    return p, a, e, i, raan, w, nu, w_true, u, lambda_true

# Vectores de ejemplo (posición y velocidad en ECI)
r_ijk = [6524.834, 6862.875, 6448.296]      # km
v_ijk = [4.901327, 5.533756, -1.976341]     # km / s

# Llamada a la función
resultados = convertir_rv_a_coe(r_ijk, v_ijk)

# Nombres de los elementos orbitales
nombres = [
    "p (semilato recto)",
    "a (semieje mayor)",
    "e (excentricidad)",
    "i (inclinación)",
    "RAAN (ascensión recta del nodo ascendente)",
    "w (argumento del perigeo)",
    "ν (anomalía verdadera)",
    "w_true (arg. perigeo verdadero)",
    "u (arg. de latitud)",
    "λ_true (longitud verdadera)"
]

# Impresión ordenada
for nombre, valor in zip(nombres, resultados):
    print(f"> {nombre:45s}: {valor}")
