import numpy as np
import matplotlib.pyplot as plt
from sgp4.api import Satrec
from sgp4.conveniences import jday
from datetime import datetime, timedelta, timezone

mu = 398600  # km^3/s^2

# Convierte epoch del TLE a datetime
def epoch_a_datetime(satelite):
    jd = satelite.jdsatepoch + satelite.jdsatepochF
    unix_time = (jd - 2440587.5) * 86400
    return datetime.fromtimestamp(unix_time, tz=timezone.utc)

# Calcula elementos orbitales clásicos
def calcular_elementos_orbitales(r, v):
    r = np.array(r, dtype=float)
    v = np.array(v, dtype=float)

    norma_r = np.linalg.norm(r)
    norma_v = np.linalg.norm(v)
    h_vec = np.cross(r, v)
    h = np.linalg.norm(h_vec)
    energia = norma_v**2 / 2 - mu / norma_r
    a = -mu / (2 * energia)
    e_vec = np.cross(v, h_vec) / mu - r / norma_r
    e = np.linalg.norm(e_vec)
    i = np.arccos(h_vec[2] / h)

    n_vec = np.cross([0, 0, 1], h_vec)
    n = np.linalg.norm(n_vec)
    RAAN = np.arccos(n_vec[0] / n) if n > 1e-8 else 0.0
    if n > 1e-8 and n_vec[1] < 0:
        RAAN = 2 * np.pi - RAAN

    w = 0.0
    if n > 1e-8 and e > 1e-8:
        w = np.arccos(np.dot(n_vec, e_vec) / (n * e))
        if e_vec[2] < 0:
            w = 2 * np.pi - w
    elif e > 1e-8:
        w = np.arctan2(e_vec[1], e_vec[0])

    f = 0.0
    if e > 1e-8:
        f = np.arccos(np.dot(e_vec, r) / (e * norma_r))
        if np.dot(r, v) < 0:
            f = 2 * np.pi - f

    return {'a': a, 'e': e, 'i': i, 'RAAN': RAAN, 'w': w, 'f': f}

# Propagación general
def propagar_desde_tle(linea1, linea2, fecha_final, paso_segundos):
    satelite = Satrec.twoline2rv(linea1, linea2)
    t_inicial = epoch_a_datetime(satelite)
    print(f"> Epoch del TLE: {t_inicial.isoformat()}")
    pasos = int((fecha_final - t_inicial).total_seconds() // paso_segundos) + 1

    resultados = []
    for paso in range(pasos):
        tiempo = t_inicial + timedelta(seconds=paso * paso_segundos)
        jd, fr = jday(tiempo.year, tiempo.month, tiempo.day,
                      tiempo.hour, tiempo.minute,
                      tiempo.second + tiempo.microsecond / 1e6)
        err, r, v = satelite.sgp4(jd, fr)
        if err != 0:
            print(f"> [!] Error {err} en t = {tiempo}")
            continue
        elementos = calcular_elementos_orbitales(r, v)
        resultados.append((tiempo, elementos))
    return resultados

# Graficar resultados
def graficar_elementos(resultados):
    tiempos = [t for t, _ in resultados]
    a = [e['a'] for _, e in resultados]
    e = [e['e'] for _, e in resultados]
    i = [np.degrees(e['i']) for _, e in resultados]
    RAAN = [np.degrees(e['RAAN']) for _, e in resultados]
    w = [np.degrees(e['w']) for _, e in resultados]
    f = [np.degrees(e['f']) for _, e in resultados]

    plt.figure(figsize=(12, 10))

    plt.subplot(3, 2, 1)
    plt.plot(tiempos, a)
    plt.title('Semieje mayor (a) [km]')
    plt.grid()

    plt.subplot(3, 2, 2)
    plt.plot(tiempos, e)
    plt.title('Excentricidad (e)')
    plt.grid()

    plt.subplot(3, 2, 3)
    plt.plot(tiempos, i)
    plt.title('Inclinación (i) [°]')
    plt.grid()

    plt.subplot(3, 2, 4)
    plt.plot(tiempos, RAAN)
    plt.title('RAAN [°]')
    plt.grid()

    plt.subplot(3, 2, 5)
    plt.plot(tiempos, w)
    plt.title('Argumento del perigeo (ω) [°]')
    plt.grid()

    plt.subplot(3, 2, 6)
    plt.plot(tiempos, f)
    plt.title('Anomalía verdadera (f) [°]')
    plt.grid()

    plt.tight_layout()
    plt.show()

# Prueba para el TLE más reciente del SAOCOM 1B de CelesTrak
if __name__ == "__main__":
    LINEA_1 = "1 46265U 20059A   25125.80864947  .00000729  00000+0  98109-4 0  9998"
    LINEA_2 = "2 46265  97.8894 312.5046 0001574  88.3358 271.8036 14.82149559253154"

    # Propagar 5 órbitas (~97.14 min x 5)
    duracion_segundos = int(97.14 * 60 * 5)
    paso_segundos = 60

    sat = Satrec.twoline2rv(LINEA_1, LINEA_2)
    fecha_final = epoch_a_datetime(sat) + timedelta(seconds=duracion_segundos)

    resultados = propagar_desde_tle(LINEA_1, LINEA_2, fecha_final, paso_segundos)
    graficar_elementos(resultados)