import numpy as np
import matplotlib.pyplot as plt

def pn0_plot_response(fsampling, Ym_pp, Ya_pp, Ym_pn, Ya_pn, Ym_np, Ya_np, Ym_nn, Ya_nn, fd0, Y_mag, Y_rad):
    # Convertir magnitudes a dB
    Ym_pp = 20 * np.log10(Ym_pp)
    Ym_pn = 20 * np.log10(Ym_pn)
    Ym_np = 20 * np.log10(Ym_np)
    Ym_nn = 20 * np.log10(Ym_nn)
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    # Ajustar los límites del eje Y para magnitud y fase
    mag_min_pp = np.nanmin(Ym_pp)
    mag_max_pp = np.nanmax(Ym_pp)
    phase_min_pp = np.nanmin(Ya_pp)
    phase_max_pp = np.nanmax(Ya_pp)
    
    mag_min_np = np.nanmin(Ym_np)
    mag_max_np = np.nanmax(Ym_np)
    phase_min_np = np.nanmin(Ya_np)
    phase_max_np = np.nanmax(Ya_np)
    
    mag_min_pn = np.nanmin(Ym_pn)
    mag_max_pn = np.nanmax(Ym_pn)
    phase_min_pn = np.nanmin(Ya_pn)
    phase_max_pn = np.nanmax(Ya_pn)
    
    mag_min_nn = np.nanmin(Ym_nn)
    mag_max_nn = np.nanmax(Ym_nn)
    phase_min_nn = np.nanmin(Ya_nn)
    phase_max_nn = np.nanmax(Ya_nn)
    
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
    axs[0, 0].semilogx(fsampling, Ym_pp, label='Theoretical response', color='k')
    axs[0, 0].set_title(r'$Y_{pp}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].set_ylim([mag_min_pp - 3, mag_max_pp + 3])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 0].semilogx(fsampling, Ya_pp, label='Theoretical response', color='k')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].set_ylim([phase_min_pp - 3, phase_max_pp + 3])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')
    
    axs[0, 1].semilogx(fsampling, Ym_pn, label='Theoretical response', color='k')
    axs[0, 1].set_title(r'$Y_{pn}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].set_ylim([mag_min_pn - 3, mag_max_pn + 3])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 1].semilogx(fsampling, Ya_pn, label='Theoretical response', color='k')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].set_ylim([phase_min_pn - 3, phase_max_pn + 3])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')
    
    axs[2, 0].semilogx(fsampling, Ym_np, label='Theoretical response', color='k')
    axs[2, 0].set_title(r'$Y_{np}(s)$')
    axs[2, 0].set_ylabel('Magnitude (dB)')
    axs[2, 0].set_xlim([x_low_axis, x_up_axis])
    axs[2, 0].set_ylim([mag_min_np - 3, mag_max_np + 3])
    axs[2, 0].grid(True, which='major', axis='both')
    axs[2, 0].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 0].semilogx(fsampling, Ya_np, label='Theoretical response', color='k')
    axs[3, 0].set_ylabel('Phase (Degrees)')
    axs[3, 0].set_xlabel('Frequency (Hz)')
    axs[3, 0].set_xlim([x_low_axis, x_up_axis])
    axs[3, 0].set_ylim([phase_min_np - 3, phase_max_np + 3])
    axs[3, 0].grid(True, which='major', axis='both')
    axs[3, 0].grid(True, which='minor', axis='x')
    
    axs[2, 1].semilogx(fsampling, Ym_nn, label='Theoretical response', color='k')
    axs[2, 1].set_title(r'$Y_{nn}(s)$')
    axs[2, 1].set_xlim([x_low_axis, x_up_axis])
    axs[2, 1].set_ylim([mag_min_nn - 3, mag_max_nn + 3])
    axs[2, 1].grid(True, which='major', axis='both')
    axs[2, 1].grid(True, which='minor', axis='x')
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[3, 1].semilogx(fsampling, Ya_nn, label='Theoretical response', color='k')
    axs[3, 1].set_xlabel('Frequency (Hz)')
    axs[3, 1].set_xlim([x_low_axis, x_up_axis])
    axs[3, 1].set_ylim([phase_min_nn - 3, phase_max_nn + 3])
    axs[3, 1].grid(True, which='major', axis='both')
    axs[3, 1].grid(True, which='minor', axis='x')

    # Agregar puntos rojos para "pn frequency scan"
    frequencies = list(Y_mag.keys())
    for freq in frequencies:
        Y_pn0_mag = Y_mag[freq]
        Y_pn0_angle = Y_rad[freq]

        axs[0, 0].semilogx(freq, 20 * np.log10(Y_pn0_mag[0, 0]), 'ro', label='pn frequency scan')
        axs[0, 1].semilogx(freq, 20 * np.log10(Y_pn0_mag[0, 1]), 'ro', label='pn frequency scan')
        axs[2, 0].semilogx(freq, 20 * np.log10(Y_pn0_mag[1, 0]), 'ro', label='pn frequency scan')
        axs[2, 1].semilogx(freq, 20 * np.log10(Y_pn0_mag[1, 1]), 'ro', label='pn frequency scan')

        axs[1, 0].semilogx(freq, (180 / np.pi) * Y_pn0_angle[0, 0], 'ro', label='pn frequency scan')
        axs[1, 1].semilogx(freq, (180 / np.pi) * Y_pn0_angle[0, 1], 'ro', label='pn frequency scan')
        axs[3, 0].semilogx(freq, (180 / np.pi) * Y_pn0_angle[1, 0], 'ro', label='pn frequency scan')
        axs[3, 1].semilogx(freq, (180 / np.pi) * Y_pn0_angle[1, 1], 'ro', label='pn frequency scan')

    axs[0, 1].legend(['Theoretical response', 'pn frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()
    