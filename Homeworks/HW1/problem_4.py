import numpy as np
import math
import matplotlib.pyplot as plt


def main():
    e_n = np.linspace(0, 1, 6)
    E_range = np.linspace(0, 2*math.pi, 100,000)

    # List of lists for M values across range of eccentricities
    M_ranges = []

    # Appending a list of M values to M_ranges for each value of eccentricity
    for i in range(6):
        M_ranges.append([E - e_n[i] * math.sin(float(E)) for E in E_range ])

    fig1, ax1 = plt.subplots()
    colours = ['r', 'b', 'c', 'g', 'm', 'k']

    # Loop zipping together e, M, and the colour for each subplot
    for colour, M_range, e in zip(colours, M_ranges, e_n):
        ax1.plot(E_range, M_range, colour, linewidth=0.5, label=f"e of {e:.1f}")

    # Labels
    ax1.set_ylabel("Mean anomaly")
    ax1.set_xlabel("Eccentric anomaly (radians)")
    ax1.set(title="Mean anomaly versus eccentric anomaly")
    ax1.legend()

    fig1.savefig("E_versus_M_plot_keplers.png")
    plt.show()
    return 0


if __name__ == '__main__':
    main()
