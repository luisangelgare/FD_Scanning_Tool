# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 13:02:44 2024

@author: Luis Angel
"""

import numpy as np
import matplotlib.pyplot as plt

#    # Función auxiliar para configurar grid y minor grid solo en el eje x
# def configure_grid_x(ax):
#     ax.grid(True, which='major', axis='both')  # Grilla mayor para ambos ejes
#     ax.xaxis.set_minor_locator(plt.AutoMinorLocator())  # Activar ticks menores solo en el eje x
#     ax.grid(which='minor', axis='x', linestyle=':', linewidth='0.5', color='gray')  # Grilla menor solo en x

def ABC_plot_response(fsampling, Fs_mag, Fs_ang, fd0, Y_mag, Y_rad):
    """
    Función para graficar la magnitud en dB y fase en grados vs frecuencia, 
    junto con un escaneo de frecuencias para el sistema ABC.

    Parámetros:
    - fsampling: Vector de frecuencias de muestreo (Hz)
    - Fs_mag: Magnitud de la función de transferencia
    - Fs_ang: Ángulo (fase) de la función de transferencia en radianes
    - fd0: Frecuencia límite inferior y superior para graficar
    - Y_mag: Diccionario con las magnitudes Y_ABC de las matrices
    - Y_rad: Diccionario con las fases Y_ABC en radianes de las matrices
    """
    
    # Convertir magnitud a decibelios y ángulo a grados
    Fs_mag_dB = 20 * np.log10(Fs_mag)
    Fs_ang_deg = (180 / np.pi) * Fs_ang
    
    # Definir límites del eje X (frecuencia)
    x_low_axis = fd0[0]
    x_up_axis = fd0[-1]
    
    # Ajustar los límites del eje Y para magnitud y fase
    mag_min = np.min(Fs_mag_dB[1:])
    mag_max = np.max(Fs_mag_dB)
    phase_min = np.min(Fs_ang_deg)
    phase_max = np.max(Fs_ang_deg)
    
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
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    # Top row - Magnitude in dB vs log(fsampling)
    axs[0, 0].semilogx(fsampling, Fs_mag_dB, label='Theoretical response', color='k')
    axs[0, 0].set_title(r'$Y_{AA}(s)$')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].set_xlim([x_low_axis, x_up_axis])
    axs[0, 0].set_ylim([mag_min - 3, mag_max + 3])
    axs[0, 0].grid(True, which='major', axis='both')
    axs[0, 0].grid(True, which='minor', axis='x')

    axs[0, 1].semilogx(fsampling, Fs_mag_dB, label='Theoretical response', color='k')
    axs[0, 1].set_title(r'$Y_{BB}(s)$')
    axs[0, 1].set_xlim([x_low_axis, x_up_axis])
    axs[0, 1].set_ylim([mag_min - 3, mag_max + 3])
    axs[0, 1].grid(True, which='major', axis='both')
    axs[0, 1].grid(True, which='minor', axis='x')

    axs[0, 2].semilogx(fsampling, Fs_mag_dB, label='Theoretical response', color='k')
    axs[0, 2].set_title(r'$Y_{CC}(s)$')
    axs[0, 2].set_xlim([x_low_axis, x_up_axis])
    axs[0, 2].set_ylim([mag_min - 3, mag_max + 3])
    axs[0, 2].grid(True, which='major', axis='both')
    axs[0, 2].grid(True, which='minor', axis='x')
    
    # Bottom row - Phase in degrees vs log(fsampling)
    axs[1, 0].semilogx(fsampling, Fs_ang_deg, label='Phase (degrees)', color='k')
    axs[1, 0].set_xlabel('Frequency (Hz)')
    axs[1, 0].set_ylabel('Phase (Degrees)')
    axs[1, 0].set_xlim([x_low_axis, x_up_axis])
    axs[1, 0].set_ylim([phase_min - 3, phase_max + 3])
    axs[1, 0].grid(True, which='major', axis='both')
    axs[1, 0].grid(True, which='minor', axis='x')

    axs[1, 1].semilogx(fsampling, Fs_ang_deg, label='Phase (degrees)', color='k')
    axs[1, 1].set_xlabel('Frequency (Hz)')
    axs[1, 1].set_xlim([x_low_axis, x_up_axis])
    axs[1, 1].set_ylim([phase_min - 3, phase_max + 3])
    axs[1, 1].grid(True, which='major', axis='both')
    axs[1, 1].grid(True, which='minor', axis='x')

    axs[1, 2].semilogx(fsampling, Fs_ang_deg, label='Phase (degrees)', color='k')
    axs[1, 2].set_xlabel('Frequency (Hz)')
    axs[1, 2].set_xlim([x_low_axis, x_up_axis])
    axs[1, 2].set_ylim([phase_min - 3, phase_max + 3])
    axs[1, 2].grid(True, which='major', axis='both')
    axs[1, 2].grid(True, which='minor', axis='x')

    # Agregar puntos rojos para "ABC frequency scan"
    frequencies = list(Y_mag.keys())
    for freq in frequencies:
        Y_ABC_mag = Y_mag[freq]
        Y_ABC_angle = Y_rad[freq]

        axs[0, 0].semilogx(freq, 20 * np.log10(Y_ABC_mag[0, 0]), 'ro', label='ABC frequency scan')
        axs[0, 1].semilogx(freq, 20 * np.log10(Y_ABC_mag[1, 1]), 'ro', label='ABC frequency scan')
        axs[0, 2].semilogx(freq, 20 * np.log10(Y_ABC_mag[2, 2]), 'ro', label=f'{freq} Hz')

        axs[1, 0].semilogx(freq, (180 / np.pi) * Y_ABC_angle[0, 0], 'ro', label='ABC frequency scan')
        axs[1, 1].semilogx(freq, (180 / np.pi) * Y_ABC_angle[1, 1], 'ro', label='ABC frequency scan')
        axs[1, 2].semilogx(freq, (180 / np.pi) * Y_ABC_angle[2, 2], 'ro', label='ABC frequency scan')

    axs[0, 2].legend(['Theoretical response', 'ABC frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()

def ABC_plot(fsampling, fd0, Y_mag, Y_rad):
    """
    Función para graficar la magnitud en dB y fase en grados vs frecuencia, 
    junto con un escaneo de frecuencias para el sistema ABC.

    Parámetros:
    - fsampling: Vector de frecuencias de muestreo (Hz)
    - fd0: Frecuencia límite inferior y superior para graficar
    - Y_mag: Diccionario con las magnitudes Y_ABC de las matrices
    - Y_rad: Diccionario con las fases Y_ABC en radianes de las matrices
    """
    
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
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    
    # Escaneo de frecuencias para ABC
    frequencies = list(Y_mag.keys())
    
    # Graficar la magnitud en dB vs log(fsampling)
    for freq in frequencies:
        Y_ABC_mag = Y_mag[freq]
        Y_ABC_angle = Y_rad[freq]

        axs[0, 0].semilogx(freq, 20 * np.log10(Y_ABC_mag[0, 0]), 'r-', label='ABC frequency scan')
        axs[0, 0].set_title(r'$Y_{AA}(s)$')
        axs[0, 0].set_ylabel('Magnitude (dB)')
        axs[0, 0].set_xlim([x_low_axis, x_up_axis])
        axs[0, 0].grid(True, which='major', axis='both')
        axs[0, 0].grid(True, which='minor', axis='x')

        axs[0, 1].semilogx(freq, 20 * np.log10(Y_ABC_mag[1, 1]), 'r-', label='ABC frequency scan')
        axs[0, 1].set_title(r'$Y_{BB}(s)$')
        axs[0, 1].set_xlim([x_low_axis, x_up_axis])
        axs[0, 1].grid(True, which='major', axis='both')
        axs[0, 1].grid(True, which='minor', axis='x')

        axs[0, 2].semilogx(freq, 20 * np.log10(Y_ABC_mag[2, 2]), 'r-', label=f'{freq} Hz')
        axs[0, 2].set_title(r'$Y_{CC}(s)$')
        axs[0, 2].set_xlim([x_low_axis, x_up_axis])
        axs[0, 2].grid(True, which='major', axis='both')
        axs[0, 2].grid(True, which='minor', axis='x')

        # Graficar la fase en grados vs log(fsampling)
        axs[1, 0].semilogx(freq, (180 / np.pi) * Y_ABC_angle[0, 0], 'r-', label='ABC frequency scan')
        axs[1, 0].set_xlabel('Frequency (Hz)')
        axs[1, 0].set_ylabel('Phase (Degrees)')
        axs[1, 0].set_xlim([x_low_axis, x_up_axis])
        axs[1, 0].grid(True, which='major', axis='both')
        axs[1, 0].grid(True, which='minor', axis='x')

        axs[1, 1].semilogx(freq, (180 / np.pi) * Y_ABC_angle[1, 1], 'r-', label='ABC frequency scan')
        axs[1, 1].set_xlabel('Frequency (Hz)')
        axs[1, 1].set_xlim([x_low_axis, x_up_axis])
        axs[1, 1].grid(True, which='major', axis='both')
        axs[1, 1].grid(True, which='minor', axis='x')

        axs[1, 2].semilogx(freq, (180 / np.pi) * Y_ABC_angle[2, 2], 'r-', label='ABC frequency scan')
        axs[1, 2].set_xlabel('Frequency (Hz)')
        axs[1, 2].set_xlim([x_low_axis, x_up_axis])
        axs[1, 2].grid(True, which='major', axis='both')
        axs[1, 2].grid(True, which='minor', axis='x')

    axs[0, 2].legend(['ABC frequency scan'], loc='lower left')

    # Ajuste final
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    plt.show()

