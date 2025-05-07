import numpy as np
import matplotlib.pyplot as plt

def ae2Eh(a, e, mu):
    if e < 0:
        raise ValueError("> La excentricidad no puede ser negativa.")

    if e == 1.0:
        # Órbita parabólica: energía = 0, a no se usa realmente
        h = np.sqrt(mu * a)  # a se interpreta como p (semi-latus rectum)
        epsilon = 0.0
        return epsilon, h

    if (e < 1 and a <= 0) or (e > 1 and a >= 0):
        raise ValueError("> La combinación de a y e no representa una órbita válida.")

    epsilon = -mu / (2 * a)
    h = np.sqrt(mu * abs(a) * abs(1 - e**2))
    return epsilon, h

def Eh2ae(epsilon, h, mu):
    if h <= 0:
        raise ValueError("> El momento angular debe ser positivo.")

    if epsilon == 0.0:
        # Órbita parabólica: a infinito, usamos p = h² / μ
        p = h**2 / mu
        a = p  # usamos p como pseudo-a
        e = 1.0
        return a, e

    if epsilon > 0:
        a = -mu / (2 * epsilon)  # para hiperbólicas: a < 0
    else:
        a = -mu / (2 * epsilon)  # elíptica: a > 0

    e2 = 1 - (h**2) / (mu * a)
    if e2 < 0:
        raise ValueError("> Los valores no representan una órbita válida.")

    e = np.sqrt(e2)
    return a, e

def plot_orbit(a, e, mu):
    epsilon, h = ae2Eh(a, e, mu)

    fig, ax = plt.subplots(figsize=(9, 6))
    theta = np.linspace(-np.pi + 0.01, np.pi - 0.01, 600)  # más amplio para ver mejor órbitas abiertas

    if e == 1.0:
        # Parabólica: usamos p = h² / μ
        p = h**2 / mu
        r = p / (1 + np.cos(theta))
        tipo = "Parabólica"
    else:
        p = a * (1 - e**2)
        r = p / (1 + e * np.cos(theta))
        tipo = "Elíptica" if e < 1 else "Hiperbólica"

    r[r > 1e8] = np.nan  # limitar extremos hiperbólicos/parabólicos

    x = r * np.cos(theta) / 1e3  # en km
    y = r * np.sin(theta) / 1e3

    ax.plot(x, y, color='purple', linewidth=2, label='Trayectoria orbital')
    ax.plot(0, 0, 'o', color='red', markersize=8, label='Cuerpo central')

    ax.set_aspect('equal')
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_title(f'>> Órbita {tipo} <<', fontsize=14)
    ax.grid(True, linestyle='--', alpha=0.5)

    texto = f"ε = {epsilon:.2e} J/kg\nh = {h:.2e} m²/s"
    ax.text(0.02, 0.98, texto, transform=ax.transAxes,
            fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.85))

    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=True)
    plt.tight_layout()
    plt.show()

def menu():
    MU = 3.986e14  # Parámetro gravitacional de la Tierra [m^3/s^2]

    while True:
        print("\n=========================== MENÚ ===========================")
        print("< 1 > Ingresar semi-eje mayor y excentricidad (a, e)")
        print("< 2 > Ingresar energía específica y momento angular (ε, h)")
        print("< 3 > Salir del programa")
        print("============================================================")

        opcion = input("\n> Seleccione una opción: ")

        try:
            if opcion == "1":
                a = float(input("> Ingrese el semi-eje mayor [m]: "))
                e = float(input("> Ingrese la excentricidad: "))
                plot_orbit(a, e, MU)

            elif opcion == "2":
                epsilon = float(input("> Ingrese la energía específica [J/kg]: "))
                h = float(input("> Ingrese el momento angular específico [m²/s]: "))
                a, e = Eh2ae(epsilon, h, MU)
                plot_orbit(a, e, MU)

            elif opcion == "3":
                print(">> ¡Gracias por utilizar el programa!")
                break

            else:
                print(">> Opción inválida. Intente nuevamente.")

        except ValueError as err:
            print(f"Error: {err}")
        except NotImplementedError as err:
            print(f"Nota: {err}")

if __name__ == "__main__":
    menu()