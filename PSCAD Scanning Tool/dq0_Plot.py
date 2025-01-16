import numpy as np
import matplotlib.pyplot as plt

def dq0_plot_response(fsampling, Ym_qq, Ya_qq, Ym_qd, Ya_qd, Ym_dq, Ya_dq, Ym_dd, Ya_dd, fd0, Y_mag, Y_rad):
    # Convertir magnitudes a dB
    Ym_qq = 20 * np.log10(Ym_qq)
    Ym_qd = 20 * np.log10(Ym_qd)
    Ym_dq = 20 * np.log10(Ym_dq)
    Ym_dd = 20 * np.log10(Ym_dd)
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    # Ajustar los límites del eje Y para magnitud y fase
    mag_min_qq = np.nanmin(Ym_qq)
    mag_max_qq = np.nanmax(Ym_qq)
    phase_min_qq = np.nanmin(Ya_qq)
    phase_max_qq = np.nanmax(Ya_qq)
    
    mag_min_dq = np.nanmin(Ym_dq)
    mag_max_dq = np.nanmax(Ym_dq)
    phase_min_dq = np.nanmin(Ya_dq)
    phase_max_dq = np.nanmax(Ya_dq)
    
    mag_min_qd = np.nanmin(Ym_qd)
    mag_max_qd = np.nanmax(Ym_qd)
    phase_min_qd = np.nanmin(Ya_qd)
    phase_max_qd = np.nanmax(Ya_qd)
    
    mag_min_dd = np.nanmin(Ym_dd)
    mag_max_dd = np.nanmax(Ym_dd)
    phase_min_dd = np.nanmin(Ya_dd)
    phase_max_dd = np.nanmax(Ya_dd)
    
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
    # Top row - Magnitude in dB vs log(fsampling)
    axs[0, 0].semilogx(fsampling, Ym_qq, label='Linear state space', color='k')
    axs[0, 0].set_title(r'$Y_{qq}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].set_ylim([mag_min_qq - 3, mag_max_qq + 3])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 0].semilogx(fsampling, Ya_qq, label='Linear state space', color='k')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].set_ylim([phase_min_qq - 3, phase_max_qq + 3])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')
    
    axs[0, 1].semilogx(fsampling, Ym_qd, label='Linear state space', color='k')
    axs[0, 1].set_title(r'$Y_{qd}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].set_ylim([mag_min_qd - 3, mag_max_qd + 3])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 1].semilogx(fsampling, Ya_qd, label='Linear state space', color='k')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].set_ylim([phase_min_qd - 3, phase_max_qd + 3])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')
    
    axs[2, 0].semilogx(fsampling, Ym_dq, label='Linear state space', color='k')
    axs[2, 0].set_title(r'$Y_{dq}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].set_ylim([mag_min_dq - 3, mag_max_dq + 3])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 0].semilogx(fsampling, Ya_dq, label='Linear state space', color='k')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].set_ylim([phase_min_dq - 3, phase_max_dq + 3])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')
    
    axs[2, 1].semilogx(fsampling, Ym_dd, label='Linear state space', color='k')
    axs[2, 1].set_title(r'$Y_{qd}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].set_ylim([mag_min_dd - 3, mag_max_dd + 3])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 1].semilogx(fsampling, Ya_dd, label='Linear state space', color='k')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].set_ylim([phase_min_dd - 3, phase_max_dd + 3])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    # Agregar puntos rojos para "dq0 frequency scan"
    frequencies = list(Y_mag.keys())
    for freq in frequencies:
        Y_qd0_mag = Y_mag[freq]
        Y_qd0_angle = Y_rad[freq]

        axs[0, 0].semilogx(freq, 20 * np.log10(Y_qd0_mag[0, 0]), 'ro', label='dq0 frequency scan')
        axs[0, 1].semilogx(freq, 20 * np.log10(Y_qd0_mag[0, 1]), 'ro', label='dq0 frequency scan')
        axs[2, 0].semilogx(freq, 20 * np.log10(Y_qd0_mag[1, 0]), 'ro', label='dq0 frequency scan')
        axs[2, 1].semilogx(freq, 20 * np.log10(Y_qd0_mag[1, 1]), 'ro', label='dq0 frequency scan')

        axs[1, 0].semilogx(freq, (180 / np.pi) * Y_qd0_angle[0, 0], 'ro', label='dq0 frequency scan')
        axs[1, 1].semilogx(freq, (180 / np.pi) * Y_qd0_angle[0, 1], 'ro', label='dq0 frequency scan')
        axs[3, 0].semilogx(freq, (180 / np.pi) * Y_qd0_angle[1, 0], 'ro', label='dq0 frequency scan')
        axs[3, 1].semilogx(freq, (180 / np.pi) * Y_qd0_angle[1, 1], 'ro', label='dq0 frequency scan')

    axs[0, 1].legend(['Linear state space', 'dq0 frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()