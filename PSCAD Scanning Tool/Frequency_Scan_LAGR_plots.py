# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 11:47:24 2024

@author: Luis Angel
"""
import matplotlib.pyplot as plt

def plot_Yabc(...):
    # Load the space-separated file into a pandas DataFrame with regex for delimiter
    df2 = pd.read_csv('Frequency_Scan_LAGR.gf81/Ys.csv', delimiter=',', header=None, engine='python')
    # Renombrar las columnas de acuerdo a su significado
    df2.columns = ['fw_str', 'Ys_str']
    # Convertir las columnas de strings complejas a números complejos
    fw = df2['fw_str'].apply(parse_complex).apply(lambda x: x.real).to_numpy()  # Extraer solo la parte real de fw
    Fs = df2['Ys_str'].apply(parse_complex).to_numpy()  # Convertir la segunda columna a números complejos
    Fs_mag=np.abs(Fs)
    Fs_ang=np.angle(Fs)
    jw1_imag = np.real(fw)/(2*np.pi)
    
    # Tiempo total de ejecución
    execution_time = end_time - start_time
    print(f"Tiempo de ejecución: {execution_time} segundos")
    
    # Convert magnitude to decibels (dB)
    Fs_mag_dB = 20 * np.log10(Fs_mag)
    Fs_ang_deg=(180/np.pi)*Fs_ang
    # Plot 6 subplots (3 on the top row, 3 on the bottom row)
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    # Top row - Magnitude in dB vs log(jw1_imag)
    axs[0, 0].plot(np.log10(jw1_imag), Fs_mag_dB, label='Magnitude (dB)', color='k')
    axs[0, 0].set_title('Y_A (s)')
    # axs[0, 0].set_xlabel('Frequency in log10(Hz)')
    axs[0, 0].set_ylabel('Magnitude (dB)')
    axs[0, 0].grid(True)
    
    axs[0, 1].plot(np.log10(jw1_imag), Fs_mag_dB, label='Magnitude (dB)', color='k')
    axs[0, 1].set_title('Y_B (s)')
    # axs[0, 1].set_xlabel('log10(Frequency)')
    # axs[0, 1].set_ylabel('Magnitude (dB)')
    axs[0, 1].grid(True)
    
    axs[0, 2].plot(np.log10(jw1_imag), Fs_mag_dB, label='Magnitude (dB)', color='k')
    axs[0, 2].set_title('Y_C (s)')
    # axs[0, 2].set_xlabel('log10(Frequency)')
    # axs[0, 2].set_ylabel('Magnitude (dB)')
    axs[0, 2].grid(True)
    
    # Bottom row - Phase in Radians vs log(jw1_imag)
    axs[1, 0].plot(np.log10(jw1_imag), Fs_ang_deg, label='Phase (degrees)', color='k')
    # axs[1, 0].set_title('Phase vs log10(jw1_imag)')
    axs[1, 0].set_xlabel('Frequency in log10(Hz)')
    axs[1, 0].set_ylabel('Phase (Radians)')
    axs[1, 0].grid(True)
    
    axs[1, 1].plot(np.log10(jw1_imag), Fs_ang_deg, label='Phase (degrees)', color='k')
    # axs[1, 1].set_title('Phase vs log10(jw1_imag)')
    axs[1, 1].set_xlabel('Frequency in log10(Hz)')
    # axs[1, 1].set_ylabel('Phase (Radians)')
    axs[1, 1].grid(True)
    
    axs[1, 2].plot(np.log10(jw1_imag), Fs_ang_deg, label='Phase (degrees)', color='k')
    # axs[1, 2].set_title('Phase vs log10(jw1_imag)')
    axs[1, 2].set_xlabel('Frequency in log10(Hz)')
    # axs[1, 2].set_ylabel('Phase (Radians)')
    axs[1, 2].grid(True)
    
    frequencies = list(Y_mag.keys())
    for freq in frequencies:
        log_freq=np.log10(freq)
        # Extract the 3x3 matrices for magnitude and phase (Y_mag, Y_rad)
        Y_ABC_mag = Y_mag[freq]
        Y_ABC_angle = Y_rad[freq]
        # Add points to the subplots (0, 0), (1, 1), and (2, 2)
        axs[0, 0].plot(log_freq, 20 * np.log10(Y_ABC_mag[0, 0]), 'ro', label=f'{freq} Hz')
        axs[0, 1].plot(log_freq, 20 * np.log10(Y_ABC_mag[1, 1]), 'ro', label=f'{freq} Hz')
        axs[0, 2].plot(log_freq, 20 * np.log10(Y_ABC_mag[2, 2]), 'ro', label=f'{freq} Hz')
    
        axs[1, 0].plot(log_freq, (180/np.pi)*Y_ABC_angle[0, 0], 'ro', label=f'{freq} Hz')
        axs[1, 1].plot(log_freq, (180/np.pi)*Y_ABC_angle[1, 1], 'ro', label=f'{freq} Hz')
        axs[1, 2].plot(log_freq, (180/np.pi)*Y_ABC_angle[2, 2], 'ro', label=f'{freq} Hz')
    
    # Adjust layout to prevent overlapping subplots
    plt.tight_layout()
    plt.show()