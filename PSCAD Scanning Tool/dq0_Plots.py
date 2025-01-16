import numpy as np
import matplotlib.pyplot as plt

def dq0_plot(fsampling, fd0, Y_mag, Y_rad):
    
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
    
    # Agregar puntos rojos para "dq0 frequency scan"
    frequencies = list(Y_mag.keys())

    # Create empty lists to store data for plotting outside the loop
    Yqq_mag, Yqq_phase, Yqd_mag, Yqd_phase = [], [], [], []
    Ydq_mag, Ydq_phase, Ydd_mag, Ydd_phase = [], [], [], []
    
    for freq in frequencies:
        Y_qd0_mag = Y_mag[freq]
        Y_qd0_angle = Y_rad[freq]

        # Collect the data
        Yqq_mag.append(20 * np.log10(Y_qd0_mag[0, 0]))
        Yqq_phase.append((180 / np.pi) * Y_qd0_angle[0, 0])
        Yqd_mag.append(20 * np.log10(Y_qd0_mag[0, 1]))
        Yqd_phase.append((180 / np.pi) * Y_qd0_angle[0, 1])
        Ydq_mag.append(20 * np.log10(Y_qd0_mag[1, 0]))
        Ydq_phase.append((180 / np.pi) * Y_qd0_angle[1, 0])
        Ydd_mag.append(20 * np.log10(Y_qd0_mag[1, 1]))
        Ydd_phase.append((180 / np.pi) * Y_qd0_angle[1, 1])

    # Now plot the full curves
    axs[0, 0].semilogx(frequencies, Yqq_mag, 'r-', label='dq0 frequency scan')
    axs[0, 0].set_title(r'$Y_{qq}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')

    axs[1, 0].semilogx(frequencies, Yqq_phase, 'r-', label='dq0 frequency scan')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')

    axs[0, 1].semilogx(frequencies, Yqd_mag, 'r-', label='dq0 frequency scan')
    axs[0, 1].set_title(r'$Y_{qd}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')

    axs[1, 1].semilogx(frequencies, Yqd_phase, 'r-', label='dq0 frequency scan')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')

    axs[2, 0].semilogx(frequencies, Ydq_mag, 'r-', label='dq0 frequency scan')
    axs[2, 0].set_title(r'$Y_{dq}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')

    axs[3, 0].semilogx(frequencies, Ydq_phase, 'r-', label='dq0 frequency scan')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')

    axs[2, 1].semilogx(frequencies, Ydd_mag, 'r-', label='dq0 frequency scan')
    axs[2, 1].set_title(r'$Y_{dd}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')

    axs[3, 1].semilogx(frequencies, Ydd_phase, 'r-', label='dq0 frequency scan')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    axs[0, 1].legend(['dq0 frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()