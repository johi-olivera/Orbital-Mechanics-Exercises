import math

# Datos
r_Tierra = 6378  # km
mu = 398600      # km^3/s^2

def calcular_semieje_mayor(r_perigeo, r_apogeo):
    # Devuelve el semieje mayor dada la distancia al perigeo y apogeo.
    return (r_perigeo + r_apogeo) / 2

def calcular_periodo_orbital(a):
    # Calcula el período orbital dado el semieje mayor "a"
    T = 2 * math.pi * math.sqrt(a**3 / mu)
    return T

def calcular_excentricidad(r_perigeo, r_apogeo):
    # Calcula la excentricidad de la órbita
    return (r_apogeo - r_perigeo) / (r_apogeo + r_perigeo)

if __name__ == "__main__":
    h_perigeo = 300     # km
    h_apogeo = 5000     # km

    r_perigeo = r_Tierra + h_perigeo
    r_apogeo = r_Tierra + h_apogeo

    a = calcular_semieje_mayor(r_perigeo, r_apogeo)
    T = calcular_periodo_orbital(a) / 60

    print(f"> El período orbital es aproximadamente {T:.2f} minutos")