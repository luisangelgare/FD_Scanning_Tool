import numpy as np
import matplotlib.pyplot as plt

def pn0_plot(fsampling, fd0, Y_mag, Y_rad):
    
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    
    # Configuración de estilos globales para las gráficas (incluyendo tamaño de textos)
    plt.rcParams.update({
        'axes.titlesize': 20,       # Tamaño de los títulos de los ejes
        'axes.labelsize': 20,       # Tamaño de las etiquetas de los ejes (X, Y)
        'xtick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje X
        'ytick.labelsize': 18,      # Tamaño del texto de las numeraciones en el eje Y
        'legend.fontsize': 20,      # Tamaño del texto de las leyendas
        'figure.titlesize': 20,     # Tamaño del título de la figura (si lo hubiera)
        'lines.linewidth': 2        # Ancho de las líneas de las curvas
    })
    
    # Crear subplots
    fig, axs = plt.subplots(4, 2, figsize=(15, 10))
    
    # Preparar listas para almacenar datos de magnitud y fase
    frequencies = sorted(Y_mag.keys())
    Ypp_mag, Ypn_mag, Ynp_mag, Ynn_mag = [], [], [], []
    Ypp_phase, Ypn_phase, Ynp_phase, Ynn_phase = [], [], [], []

    # Acumular valores para cada frecuencia
    for freq in frequencies:
        Y_pn0_mag = Y_mag[freq]
        Y_pn0_angle = Y_rad[freq]
        
        # Magnitudes
        Ypp_mag.append(20 * np.log10(Y_pn0_mag[0, 0]) if Y_pn0_mag[0, 0] > 0 else None)
        Ypn_mag.append(20 * np.log10(Y_pn0_mag[0, 1]) if Y_pn0_mag[0, 1] > 0 else None)
        Ynp_mag.append(20 * np.log10(Y_pn0_mag[1, 0]) if Y_pn0_mag[1, 0] > 0 else None)
        Ynn_mag.append(20 * np.log10(Y_pn0_mag[1, 1]) if Y_pn0_mag[1, 1] > 0 else None)

        # Fases
        Ypp_phase.append((180 / np.pi) * Y_pn0_angle[0, 0])
        Ypn_phase.append((180 / np.pi) * Y_pn0_angle[0, 1])
        Ynp_phase.append((180 / np.pi) * Y_pn0_angle[1, 0])
        Ynn_phase.append((180 / np.pi) * Y_pn0_angle[1, 1])

    # Graficar las magnitudes y fases acumuladas
    axs[0, 0].semilogx(frequencies, Ypp_mag, 'r-', label='pn frequency scan')
    axs[0, 0].set_title(r'$Y_{pp}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')

    axs[1, 0].semilogx(frequencies, Ypp_phase, 'r-', label='pn frequency scan')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')

    axs[0, 1].semilogx(frequencies, Ypn_mag, 'r-', label='pn frequency scan')
    axs[0, 1].set_title(r'$Y_{pn}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')

    axs[1, 1].semilogx(frequencies, Ypn_phase, 'r-', label='pn frequency scan')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')

    axs[2, 0].semilogx(frequencies, Ynp_mag, 'r-', label='pn frequency scan')
    axs[2, 0].set_title(r'$Y_{np}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')

    axs[3, 0].semilogx(frequencies, Ynp_phase, 'r-', label='pn frequency scan')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')

    axs[2, 1].semilogx(frequencies, Ynn_mag, 'r-', label='pn frequency scan')
    axs[2, 1].set_title(r'$Y_{nn}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')

    axs[3, 1].semilogx(frequencies, Ynn_phase, 'r-', label='pn frequency scan')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    axs[0, 1].legend(['pn frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()
