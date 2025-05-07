import numpy as np
import matplotlib.pyplot as plt

def ae2Eh(a, e, mu):
    """
    Calcula la energía específica (ε) y el momento angular específico (h)
    a partir del semi-eje mayor (a) y la excentricidad (e) de una órbita elíptica.

    Parámetros:
        a (float): semi-eje mayor [m]. Debe ser mayor a 0.
        e (float): excentricidad [-]. Debe estar en el rango [0, 1).
        mu (float): parámetro gravitacional [m^3/s^2].

    Retorna:
        epsilon (float): energía específica [J/kg].
        h (float): módulo del momento angular específico [m²/s].

    Restricciones:
        - Si a <= 0: no se puede definir una órbita.
        - Si e >= 1 o e < 0: no corresponde a una órbita elíptica válida.
    
    Notas:
        - Para e = 0 la órbita es circular (válida).
        - No se consideran órbitas parabólicas o hiperbólicas.
    """
    if a <= 0:
        raise ValueError("> El semi-eje mayor debe ser positivo.")
    if not (0 <= e < 1):
        raise ValueError("> La excentricidad debe estar entre 0 (incluido) y 1 (excluido).")
    
    epsilon = -mu / (2 * a)
    h = np.sqrt(mu * a * (1 - e**2))
    return epsilon, h

def Eh2ae(epsilon, h, mu):
    """
    Reconstruye el semi-eje mayor (a) y la excentricidad (e)
    a partir de la energía específica (ε) y el momento angular (h).

    Parámetros:
        epsilon (float): energía específica [J/kg]. Debe ser negativa.
        h (float): módulo del momento angular específico [m²/s]. Debe ser positivo.
        mu (float): parámetro gravitacional [m^3/s^2].

    Retorna:
        a (float): semi-eje mayor [m].
        e (float): excentricidad [-].

    Restricciones:
        - Si ε >= 0: no representa una órbita elíptica (podría ser parabólica o escape).
        - Si h <= 0: momento angular no válido.
        - Si el cálculo de e² da negativo: los valores no representan una órbita física válida.
    """
    if epsilon >= 0:
        raise ValueError("> La energía debe ser negativa para órbitas elípticas.")
    if h <= 0:
        raise ValueError("> El momento angular debe ser positivo.")
    
    a = -mu / (2 * epsilon)
    e2 = 1 - (h**2) / (mu * a)
    
    if e2 < 0:
        raise ValueError("> Los valores no representan una órbita elíptica válida.")
    
    e = np.sqrt(e2)
    return a, e

def plot_orbit(a, e, mu):
    """
    Dibuja una órbita elíptica con la trayectoria en violeta y el cuerpo central en rojo.
    Muestra una leyenda con información física y etiquetas al costado.

    Parámetros:
        a (float): semi-eje mayor [m].
        e (float): excentricidad [-].
        mu (float): parámetro gravitacional [m^3/s^2].
    """
    epsilon, h = ae2Eh(a, e, mu)
    
    # Parámetro orbital
    p = a * (1 - e**2)
    theta = np.linspace(0, 2 * np.pi, 600)
    r = p / (1 + e * np.cos(theta))
    
    # Coordenadas cartesianas (en km)
    x = r * np.cos(theta) / 1e3
    y = r * np.sin(theta) / 1e3
    
    # Gráfico
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(x, y, color='purple', linewidth=2, label='Trayectoria orbital')
    ax.plot(0, 0, 'o', color='red', markersize=8, label='Cuerpo central')
    
    ax.set_aspect('equal')
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_title('>> Órbita Elíptica <<', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.5)

    # Información de energía y momento angular
    texto = f"ε = {epsilon:.2e} J/kg\nh = {h:.2e} m²/s"
    ax.text(0.02, 0.98, texto, transform=ax.transAxes,
            fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.85))

    # Leyenda al costado
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=True)
    
    plt.tight_layout()
    plt.show()

# ==========================
#    >> Caso de prueba <<
# ==========================

if __name__ == "__main__":
    MU = 3.986e14       # Parámetro gravitacional estándar [m^3/s^2]
    a_demo = 15000e3    # 15,000 km en metros
    e_demo = 0.6        # excentricidad más alta

    ε, h = ae2Eh(a_demo, e_demo, MU)
    print(f">> [ae2Eh] a = {a_demo/1e3:.0f} km, e = {e_demo} -> ε = {ε:.2e} J/kg, h = {h:.2e} m²/s")

    a_rec, e_rec = Eh2ae(ε, h, MU)
    print(f">> [Eh2ae] ε = {ε:.2e}, h = {h:.2e} -> a = {a_rec/1e3:.0f} km, e = {e_rec:.4f}")

    plot_orbit(a_demo, e_demo, MU)