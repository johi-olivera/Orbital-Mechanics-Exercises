import math
from Ejercicio13 import (
    calcular_semieje_mayor,
    calcular_periodo_orbital,
    calcular_excentricidad,
    r_Tierra
)

# Datos
h_perigeo = 300     # km
h_apogeo = 5000     # km

r_perigeo = r_Tierra + h_perigeo
r_apogeo = r_Tierra + h_apogeo

# Cálculos orbitales
a = calcular_semieje_mayor(r_perigeo, r_apogeo)
T = calcular_periodo_orbital(a)
e = calcular_excentricidad(r_perigeo, r_apogeo)

# Tiempo transcurrido después del apogeo
delta_t = 80 * 60  # Segundos

# Tiempo en el apogeo (E = π)
E_apogeo = math.pi
t_apogeo = (T / (2 * math.pi)) * (E_apogeo - e * math.sin(E_apogeo))
t_total = t_apogeo + delta_t

# Anomalía media
n = 2 * math.pi / T  # rad/s
M = n * t_total

# Resolver ecuación de Kepler: M = E - e * sin(E)
def resolver_Kepler(tolerancia_grados):
    tolerancia_rad = math.radians(tolerancia_grados)
    E = M if e < 0.8 else math.pi
    while True:
        funcion = E - e * math.sin(E) - M
        derivada = 1 - e * math.cos(E)
        E_nueva = E - funcion / derivada
        if abs(E_nueva - E) < tolerancia_rad:
            break
        E = E_nueva
    return E

# Convertir E a ν (anomalía verdadera)
def convertir_E_a_nu(E):
    cos_nu = (math.cos(E) - e) / (1 - e * math.cos(E))
    sin_nu = (math.sqrt(1 - e**2) * math.sin(E)) / (1 - e * math.cos(E))
    nu = math.atan2(sin_nu, cos_nu)
    if nu < 0:
        nu += 2 * math.pi
    return nu

# Cálculo final
tolerancia_grados = 0.1    # Grados

E_final = resolver_Kepler(tolerancia_grados)
nu_final = convertir_E_a_nu(E_final)

nu_grados = math.degrees(nu_final)

print(f"> La anomalía verdadera 80 minutos después del apogeo es aproximadamente {nu_grados:.2f}°")