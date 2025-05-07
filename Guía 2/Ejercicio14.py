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

# Conversión de anomalía verdadera a excéntrica
def true_to_eccentric_anomaly(nu, e):
    tan_half_E = math.sqrt((1 - e) / (1 + e)) * math.tan(nu / 2)
    E = 2 * math.atan(tan_half_E)
    if E < 0:
        E += 2 * math.pi
    return E

# Anomalías verdaderas
nu1 = math.radians(90)
nu2 = math.radians(180)

# Anomalías excéctricas
E1 = true_to_eccentric_anomaly(nu1, e)
E2 = true_to_eccentric_anomaly(nu2, e)

# Tiempos desde el perigeo
t1 = (T / (2*math.pi)) * (E1 - e * math.sin(E1))
t2 = (T / (2*math.pi)) * (E2 - e * math.sin(E2))

delta_t_min = (t2 - t1) / 60    # En minutos

print(f"> El período orbital es aproximadamente {T / 60:.2f} minutos")
print(f"> El tiempo desde el semi-lactus rectum (90°) al apogeo es aproximadamente {delta_t_min:.2f} minutos")