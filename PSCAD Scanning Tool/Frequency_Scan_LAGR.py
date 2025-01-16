# -*- coding: utf-8 -*-

# % Technical University of Catalonia (UPC)
# % Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
# % Centre of Technological Innovation in Static Converters and Drives (CITCEA)
# % Doctoral Program in Electrical Engineering
# % Developed by: Luis Angel Garcia Reyes, MSc

# Frequency domain scanner

import mhi.pscad
import numpy as np
import pandas as pd
import sys
from decimal import Decimal
import time
import matplotlib.pyplot as plt
from ABC_Plot import ABC_plot_response, ABC_plot
from dq0_Plot import dq0_plot_response
from pn0_Plot import pn0_plot_response
from pn0_Plots import pn0_plot
from dq0_Plots import dq0_plot
from scipy.io import savemat


# Cartesian to polar function
def cartz2pol(z):
    r = np.abs(z)
    theta = np.angle(z)
    return r, theta

# Función para convertir cadena compleja a número complejo en Python
def parse_complex(complex_str):
    # Reemplazamos 'i' por 'j' para que Python lo interprete como complejo
    complex_str = complex_str.replace('i', 'j')
    return complex(complex_str)

# Start point
start_time = time.time()

# Open and launch PSCAD
pscad = mhi.pscad.launch()

# Find and load the PSCAD project
pscad.load(r'C:\Users\Luis Angel\Documents\Doctorado\PSCAD\Frequency Scan\Frequency_Scan_LAGR.pscx')
project_name = 'Frequency_Scan_LAGR'
project = pscad.project(project_name)
# Modify and access to the simulation parameters of PSCAD
scanner_selection=2
Tobs = 3.0   # Total simulation time 
step_time = 25 # step time for PSCAD (coeficient only)
plot_step = step_time*1 # Channel Plot Step (how many steps time are registered)
project.parameters({
    'time_duration': Tobs,
    'time_step': step_time,
    'sample_step': plot_step
})
# Get the "Main" canvas (layer or space) from PSCAD
main = project.canvas('Main')

# Setting the simulation parameters
coef_decimal = Decimal('1E-6')
delta_t = float(coef_decimal * Decimal(step_time)) # Step time in real number
samples = int(Tobs / delta_t)  # Number of samples
Delta_f0 = 1 / Tobs  # Frequency step (rad)
f1 = np.arange(0, samples) * Delta_f0  # Frequency vector
jw1 = 1j * 2 * np.pi * f1  # Complex frequency vector
t_w1 = 1.8  # Starting time window
t_w2 = 2.8  # Ending time window
Tobs0 = t_w2 - t_w1  # Total time window
fs = 1 / Tobs0  # Sampling freq1ºuency for FFT

# Power system parameters
Vbase=66E3
f0 = 50  # Base frequency
phi = (2 / 3) * np.pi  # Phase of 120º between phases
dist_act = 0.0 # Disturbance time
Vdist=0.03*Vbase

#### FREQUENCY SCAN SETTINGS #####
# Generate logarithmic spaced frequencies values 
fd0 = np.unique(np.round(np.logspace(-1, 3, 150))) 
# fd0=np.unique(np.round(np.logspace(2, 2.5, 50)))
# Find the main components by the ID in PSCAD
Scanner = project.component(502873399)
Scanner.parameters(Vbase=Vbase, fbase=f0, Rsource=1E-3, dist_time=dist_act)
# Define the dictionaries spaces for Y calculations
Y = {}
Y_mag = {}
Y_rad = {}

#-------------- ABC scanner ---------------------------------------------------

if scanner_selection ==1:
    ABC_Scan_Layer = project.set_layer('ABC_Scanner', 'enabled')
    dq0_Scan_Layer = project.set_layer('dq0_Scanner', 'disabled')
    pn0_Scan_Layer = project.set_layer('pn0_Scanner', 'disabled')
    for i in range(len(fd0)):
        # Assign the value to the variable in PSCAD
        Scanner.parameters(Vdist=Vdist, fd=fd0[i])
        # print(Scanner.parameters())
        # Define the dictionaries spaces for Y calculations)
        # Run the simulation
        project.run()
        # Load the space-separated file into a pandas DataFrame with regex for delimiter
        df = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_01.out', delimiter='\s+', header=None, engine='python')
        # Drop any columns that contain only NaN values
        df_cleaned = df.dropna(axis=1, how='all')
        # Set column names to 'time', 'A', 'B', 'C'
        df_cleaned.columns = ['time','ia', 'ib', 'ic','va', 'vb', 'vc']
        df_cleaned = df_cleaned.set_index('time')
        time_w = df_cleaned.loc[t_w1:t_w2]
        # Extract the V and I windowed signals
        va = time_w['va'].to_numpy()
        vb = time_w['vb'].to_numpy()
        vc = time_w['vc'].to_numpy()
        ia = time_w['ia'].to_numpy()
        ib = time_w['ib'].to_numpy()
        ic = time_w['ic'].to_numpy()
        # Apply the FFT to each signal
        Va = np.fft.fft(va) / len(va)
        Vb = np.fft.fft(vb) / len(vb)
        Vc = np.fft.fft(vc) / len(vc)
        Ia = np.fft.fft(ia) / len(ia)
        Ib = np.fft.fft(ib) / len(ib)
        Ic = np.fft.fft(ic) / len(ic)
        # Generating the corresponding frequencies based on the window time
        frequencies = np.fft.fftfreq(len(va), d=Tobs0)
        # Apply the unilateral sprectrum 
        positive_freqs = frequencies[frequencies >= 0]
        Va = Va[:len(positive_freqs)]
        Vb = Vb[:len(positive_freqs)]
        Vc = Vc[:len(positive_freqs)]
        Ia = Ia[:len(positive_freqs)]
        Ib = Ib[:len(positive_freqs)]
        Ic = Ic[:len(positive_freqs)]
        # Generating the frequency response at f_disturbance and build the Y matrix
        wd = round(fd0[i]/fs)
        Yaa=Ia[wd]/Va[wd]
        Yab=Ia[wd]/Vb[wd]
        Yac=Ia[wd]/Vc[wd]
        Yba=Ib[wd]/Va[wd]
        Ybb=Ib[wd]/Vb[wd]
        Ybc=Ib[wd]/Vc[wd]
        Yca=Ic[wd]/Va[wd]
        Ycb=Ic[wd]/Vb[wd]
        Ycc=Ic[wd]/Vc[wd]
        Yabc = np.array([[Yaa, Yab, Yac],
                        [Yba, Ybb, Ybc],
                        [Yca, Ycb, Ycc]])
        Y[fd0[i]]=Yabc
        [Y_mag[fd0[i]], Y_rad[fd0[i]]]=cartz2pol(Yabc)
    # Close PSCAD once the scan has finished
    pscad.quit() 
    # Load the RLC results:R = 1, L = 0.1, C = 10E-6
    # Load the space-separated file into a pandas DataFrame with regex for delimiter
    df2 = pd.read_csv('Frequency_Scan_LAGR.gf81/Ys.csv', delimiter=',', header=None, engine='python')
    # Renombrar las columnas de acuerdo a su significado
    df2.columns = ['jw_str', 'Ys_str']
    # Convertir las columnas de strings complejas a números complejos
    jw_sampling = df2['jw_str'].apply(parse_complex).apply(lambda x: x.real).to_numpy()  # Extraer solo la parte real de fw
    Fs = df2['Ys_str'].apply(parse_complex).to_numpy()  # Convertir la segunda columna a números complejos
    Fs_mag=np.abs(Fs)
    Fs_ang=np.angle(Fs)
    fsampling = np.real(jw_sampling)/(2*np.pi)
    ABC_plot_response(fsampling, Fs_mag, Fs_ang, fd0, Y_mag, Y_rad) 
     
#-------------- dq0 scanner ---------------------------------------------------

if scanner_selection ==2:
    ABC_Scan_Layer = project.set_layer('ABC_Scanner', 'disabled')
    dq0_Scan_Layer = project.set_layer('dq0_Scanner', 'enabled')
    pn0_Scan_Layer = project.set_layer('pn0_Scanner', 'disabled')
    dq0_block = project.component(1391233784)
    for i in range(len(fd0)):
# Injection in q-sequence
        # Assign the value to the variable in PSCAD
        Scanner.parameters(Vdist=Vdist, fd=fd0[i])
        dq0_block.parameters(Vdist_q=Vdist, Vdist_d=0)
        # Define the dictionaries spaces for Y calculations)
        # Run the simulation
        project.run()
        # Load the space-separated file into a pandas DataFrame with regex for delimiter
        df = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_01.out', delimiter='\s+', header=None, engine='python')
        # Drop any columns that contain only NaN values
        df_cleaned = df.dropna(axis=1, how='all')
        # Set column names to 'time', 'A', 'B', 'C'
        df_cleaned.columns = ['time','i_q', 'i_d', 'v_q','v_d']
        df_cleaned = df_cleaned.set_index('time')
        time_w = df_cleaned.loc[t_w1:t_w2]
        # Extract the V and I windowed signals
        v_q = time_w['v_q'].to_numpy()
        v_d = time_w['v_d'].to_numpy()
        i_q = time_w['i_q'].to_numpy()
        i_d = time_w['i_d'].to_numpy()
        # Apply the FFT to each signal
        Vq = np.fft.fft(v_q) / len(v_q)
        Vd = np.fft.fft(v_d) / len(v_d)
        Iq = np.fft.fft(i_q) / len(i_q)
        Id = np.fft.fft(i_d) / len(i_d)
        # Generating the corresponding frequencies based on the window time
        frequencies = np.fft.fftfreq(len(v_q), d=Tobs0)
        # Apply the unilateral sprectrum 
        positive_freqs = frequencies[frequencies >= 0]
        Vq = Vq[:len(positive_freqs)]
        Vd = Vd[:len(positive_freqs)]
        Iq = Iq[:len(positive_freqs)]
        Id = Id[:len(positive_freqs)]
        # Generating the frequency response at f_disturbance and build the Y matrix
        wd = round(fd0[i]/fs)
        Yqq=Iq[wd]/Vq[wd]
        Ydq=Id[wd]/Vq[wd]
# Injection in d-sequence
        # Assign the value to the variable in PSCAD
        dq0_block.parameters(Vdist_q=0, Vdist_d=Vdist)
        # Define the dictionaries spaces for Y calculations)
        # Run the simulation
        project.run()
        # Load the space-separated file into a pandas DataFrame with regex for delimiter
        df = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_01.out', delimiter='\s+', header=None, engine='python')
        # Drop any columns that contain only NaN values
        df_cleaned = df.dropna(axis=1, how='all')
        # Set column names to 'time', 'A', 'B', 'C'
        df_cleaned.columns = ['time','i_q', 'i_d', 'v_q','v_d']
        df_cleaned = df_cleaned.set_index('time')
        time_w = df_cleaned.loc[t_w1:t_w2]
        # Extract the V and I windowed signals
        v_q = time_w['v_q'].to_numpy()
        v_d = time_w['v_d'].to_numpy()
        i_q = time_w['i_q'].to_numpy()
        i_d = time_w['i_d'].to_numpy()
        # Apply the FFT to each signal
        Vq = np.fft.fft(v_q) / len(v_q)
        Vd = np.fft.fft(v_d) / len(v_d)
        Iq = np.fft.fft(i_q) / len(i_q)
        Id = np.fft.fft(i_d) / len(i_d)
        # Apply the unilateral sprectrum 
        Vq = Vq[:len(positive_freqs)]
        Vd = Vd[:len(positive_freqs)]
        Iq = Iq[:len(positive_freqs)]
        Id = Id[:len(positive_freqs)]
        # Generating the frequency response at f_disturbance and build the Y matrix
        Yqd=Iq[wd]/Vd[wd]
        Ydd=Id[wd]/Vd[wd]
        Yqd0 = np.array([[Yqq, Yqd],
                        [Ydq, Ydd]])
        Y[fd0[i]]=Yqd0
        [Y_mag[fd0[i]], Y_rad[fd0[i]]]=cartz2pol(Yqd0)
    # Open loop for varying the parameter value in the scanner
    # Close PSCAD once the scan has finished
    pscad.quit()    
    # Load the space-separated file into a pandas DataFrame with regex for delimiter
    df2 = pd.read_csv('Frequency_Scan_LAGR.gf81/Ymatrix.csv', delimiter=',', header=None, engine='python')
    # Renombrar las columnas de acuerdo a su significado
    df2.columns = ['jw_str', 'Ym_qq', 'Ya_qq', 'Ym_qd', 'Ya_qd', 'Ym_dq', 'Ya_dq', 'Ym_dd', 'Ya_dd',]
    # Convertir las columnas de strings complejas a números complejos
    jw_sampling = df2['jw_str'].apply(lambda x: x.real).to_numpy()  # Extraer solo la parte real de fw
    Ym_qq = df2['Ym_qq'].to_numpy()
    Ya_qq = df2['Ya_qq'].to_numpy()
    Ym_qd = df2['Ym_qd'].to_numpy()
    Ya_qd = df2['Ya_qd'].to_numpy()
    Ym_dq = df2['Ym_dq'].to_numpy()
    Ya_dq = df2['Ya_dq'].to_numpy()
    Ym_dd = df2['Ym_dd'].to_numpy()
    Ya_dd = df2['Ya_dd'].to_numpy()
    fsampling = np.real(jw_sampling)/(2*np.pi)
    # dq0_plot_response(fsampling, Ym_qq, Ya_qq, Ym_qd, Ya_qd, Ym_dq, Ya_dq, Ym_dd, Ya_dd, fd0, Y_mag, Y_rad)
    dq0_plot(fsampling, fd0, Y_mag, Y_rad)
    # Guarda los diccionarios Y_mag y Y_rad en un archivo .mat
    # output_data = {'Y_mag': Y_mag, 'Y_rad': Y_rad}
    # Especifica el nombre del archivo .mat
    # savemat('Y_matrices_bergeron.mat', output_data)
#-------------- pn0 scanner ---------------------------------------------------

if scanner_selection ==3:
    ABC_Scan_Layer = project.set_layer('ABC_Scanner', 'disabled')
    dq0_Scan_Layer = project.set_layer('dq0_Scanner', 'disabled')
    pn0_Scan_Layer = project.set_layer('pn0_Scanner', 'enabled')
    pn0_block = project.component(1518213203)
    for i in range(len(fd0)):
# Injection in p-sequence
        # Assign the value to the variable in PSCAD
        Scanner.parameters(Vdist=Vdist, fd=fd0[i])
        pn0_block.parameters(Vdist_p=Vdist, Vdist_n=0.0, Vdist_0=0.0)
        # Define the dictionaries spaces for Y calculations)
        # Run the simulation
        project.run()
        # Load the space-separated file into a pandas DataFrame with regex for delimiter
        df = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_01.out', delimiter='\s+', header=None, engine='python')
        df2 = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_02.out', delimiter='\s+', header=None, engine='python')
        # Drop any columns that contain only NaN values
        df_cleaned = df.dropna(axis=1, how='all')
        df_cleaned2 = df2.dropna(axis=1, how='all')
        # Set column names to 'time', 'Re{i0}', 'I{i0}', 'Re{ip}', 'I{ip}', 'Re{in}', 'I{in}', ...
        df_cleaned.columns = ['time','R_i_0', 'R_i_p',  'R_i_n', 'I_i_0', 'I_i_p', 'I_i_n', 'R_v_0', 'R_v_p', 'R_v_n', 'I_v_0']
        df_cleaned2.columns = ['time', 'I_v_p', 'I_v_n']
        df_cleaned = df_cleaned.set_index('time')
        df_cleaned2 = df_cleaned2.set_index('time')
        # Combinar las partes real e imaginaria en números complejos
        df_cleaned['i_0'] = df_cleaned['R_i_0'] + 1j * df_cleaned['I_i_0']
        df_cleaned['i_p'] = df_cleaned['R_i_p'] + 1j * df_cleaned['I_i_p']
        df_cleaned['i_n'] = df_cleaned['R_i_n'] + 1j * df_cleaned['I_i_n']
        df_cleaned['v_0'] = df_cleaned['R_v_0'] + 1j * df_cleaned['I_v_0']
        df_cleaned['v_p'] = df_cleaned['R_v_p'] + 1j * df_cleaned2['I_v_p']
        df_cleaned2['v_n'] = df_cleaned['R_v_n'] + 1j * df_cleaned2['I_v_n']
        # Concatenar df_cleaned2 en df_cleaned para tener todos los datos en un solo DataFrame
        df_combined = pd.concat([df_cleaned, df_cleaned2], axis=1)
        time_w = df_combined.loc[t_w1:t_w2]
        # Extract the V and I windowed signals
        v_p = time_w['v_p'].to_numpy()
        v_n = time_w['v_n'].to_numpy()
        v_0 = time_w['v_0'].to_numpy()
        i_p = time_w['i_p'].to_numpy()
        i_n = time_w['i_n'].to_numpy()
        i_0 = time_w['i_0'].to_numpy()
        # Apply the FFT to each signal
        Vp = np.fft.fft(v_p) / len(v_p)
        Vn = np.fft.fft(v_n) / len(v_n)
        V0 = np.fft.fft(v_0) / len(v_0)
        Ip = np.fft.fft(i_p) / len(i_p)
        In = np.fft.fft(i_n) / len(i_n)
        I0 = np.fft.fft(i_0) / len(i_0)
        # Generating the corresponding frequencies based on the window time
        frequencies = np.fft.fftfreq(len(v_p), d=Tobs0)
        # Apply the unilateral sprectrum 
        positive_freqs = frequencies[frequencies >= 0]
        Vp = Vp[:len(positive_freqs)]
        Vn = Vn[:len(positive_freqs)]
        V0 = V0[:len(positive_freqs)]
        Ip = Ip[:len(positive_freqs)]
        In = In[:len(positive_freqs)]
        I0 = I0[:len(positive_freqs)]
        # Generating the frequency response at f_disturbance and build the Y matrix
        wd = round(fd0[i]/fs)
        Ypp=Ip[wd]/Vp[wd]
        Ynp=In[wd]/Vp[wd]
        Y0p=I0[wd]/Vp[wd]
# Injection in n-sequence
        # Assign the value to the variable in PSCAD
        Scanner.parameters(Vdist=Vdist, fd=-1*fd0[i])
        pn0_block.parameters(Vdist_p=0, Vdist_n=Vdist, Vdist_0=0)
        # Define the dictionaries spaces for Y calculations)
        # Run the simulation
        project.run()
        # Load the space-separated file into a pandas DataFrame with regex for delimiter
        df = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_01.out', delimiter='\s+', header=None, engine='python')
        df2 = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_02.out', delimiter='\s+', header=None, engine='python')
        # Drop any columns that contain only NaN values
        df_cleaned = df.dropna(axis=1, how='all')
        df_cleaned2 = df2.dropna(axis=1, how='all')
        # Set column names to 'time', 'Re{i0}', 'I{i0}', 'Re{ip}', 'I{ip}', 'Re{in}', 'I{in}', ...
        df_cleaned.columns = ['time','R_i_0', 'R_i_p',  'R_i_n', 'I_i_0', 'I_i_p', 'I_i_n', 'R_v_0', 'R_v_p', 'R_v_n', 'I_v_0']
        df_cleaned2.columns = ['time', 'I_v_p', 'I_v_n']
        df_cleaned = df_cleaned.set_index('time')
        df_cleaned2 = df_cleaned2.set_index('time')
        # Combinar las partes real e imaginaria en números complejos
        df_cleaned['i_0'] = df_cleaned['R_i_0'] + 1j * df_cleaned['I_i_0']
        df_cleaned['i_p'] = df_cleaned['R_i_p'] + 1j * df_cleaned['I_i_p']
        df_cleaned['i_n'] = df_cleaned['R_i_n'] + 1j * df_cleaned['I_i_n']
        df_cleaned['v_0'] = df_cleaned['R_v_0'] + 1j * df_cleaned['I_v_0']
        df_cleaned['v_p'] = df_cleaned['R_v_p'] + 1j * df_cleaned2['I_v_p']
        df_cleaned2['v_n'] = df_cleaned['R_v_n'] + 1j * df_cleaned2['I_v_n']
        # Concatenar df_cleaned2 en df_cleaned para tener todos los datos en un solo DataFrame
        df_combined = pd.concat([df_cleaned, df_cleaned2], axis=1)
        time_w = df_combined.loc[t_w1:t_w2]
        # Extract the V and I windowed signals
        v_p = time_w['v_p'].to_numpy()
        v_n = time_w['v_n'].to_numpy()
        v_0 = time_w['v_0'].to_numpy()
        i_p = time_w['i_p'].to_numpy()
        i_n = time_w['i_n'].to_numpy()
        i_0 = time_w['i_0'].to_numpy()
        # Apply the FFT to each signal
        Vp = np.fft.fft(v_p) / len(v_p)
        Vn = np.fft.fft(v_n) / len(v_n)
        V0 = np.fft.fft(v_0) / len(v_0)
        Ip = np.fft.fft(i_p) / len(i_p)
        In = np.fft.fft(i_n) / len(i_n)
        I0 = np.fft.fft(i_0) / len(i_0)
        Vp = Vp[:len(positive_freqs)]
        Vn = Vn[:len(positive_freqs)]
        V0 = V0[:len(positive_freqs)]
        Ip = Ip[:len(positive_freqs)]
        In = In[:len(positive_freqs)]
        I0 = I0[:len(positive_freqs)]
        Ypn=Ip[wd]/Vn[wd]
        Ynn=In[wd]/Vn[wd]
        Y0n=I0[wd]/Vn[wd]
# Injection in 0-sequence
        # Assign the value to the variable in PSCAD
        Scanner.parameters(Vdist=Vdist, fd=fd0[i])
        pn0_block.parameters(Vdist_p=0, Vdist_n=0, Vdist_0=Vdist)
        # Define the dictionaries spaces for Y calculations)
        # Run the simulation
        project.run()
        # Load the space-separated file into a pandas DataFrame with regex for delimiter
        df = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_01.out', delimiter='\s+', header=None, engine='python')
        df2 = pd.read_csv('Frequency_Scan_LAGR.gf81/VI_abc_02.out', delimiter='\s+', header=None, engine='python')
        # Drop any columns that contain only NaN values
        df_cleaned = df.dropna(axis=1, how='all')
        df_cleaned2 = df2.dropna(axis=1, how='all')
        # Set column names to 'time', 'Re{i0}', 'I{i0}', 'Re{ip}', 'I{ip}', 'Re{in}', 'I{in}', ...
        df_cleaned.columns = ['time','R_i_0', 'R_i_p',  'R_i_n', 'I_i_0', 'I_i_p', 'I_i_n', 'R_v_0', 'R_v_p', 'R_v_n', 'I_v_0']
        df_cleaned2.columns = ['time', 'I_v_p', 'I_v_n']
        df_cleaned = df_cleaned.set_index('time')
        df_cleaned2 = df_cleaned2.set_index('time')
        # Combinar las partes real e imaginaria en números complejos
        df_cleaned['i_0'] = df_cleaned['R_i_0'] + 1j * df_cleaned['I_i_0']
        df_cleaned['i_p'] = df_cleaned['R_i_p'] + 1j * df_cleaned['I_i_p']
        df_cleaned['i_n'] = df_cleaned['R_i_n'] + 1j * df_cleaned['I_i_n']
        df_cleaned['v_0'] = df_cleaned['R_v_0'] + 1j * df_cleaned['I_v_0']
        df_cleaned['v_p'] = df_cleaned['R_v_p'] + 1j * df_cleaned2['I_v_p']
        df_cleaned2['v_n'] = df_cleaned['R_v_n'] + 1j * df_cleaned2['I_v_n']
        # Concatenar df_cleaned2 en df_cleaned para tener todos los datos en un solo DataFrame
        df_combined = pd.concat([df_cleaned, df_cleaned2], axis=1)
        time_w = df_combined.loc[t_w1:t_w2]
        # Extract the V and I windowed signals
        v_p = time_w['v_p'].to_numpy()
        v_n = time_w['v_n'].to_numpy()
        v_0 = time_w['v_0'].to_numpy()
        i_p = time_w['i_p'].to_numpy()
        i_n = time_w['i_n'].to_numpy()
        i_0 = time_w['i_0'].to_numpy()
        # Apply the FFT to each signal
        Vp = np.fft.fft(v_p) / len(v_p)
        Vn = np.fft.fft(v_n) / len(v_n)
        V0 = np.fft.fft(v_0) / len(v_0)
        Ip = np.fft.fft(i_p) / len(i_p)
        In = np.fft.fft(i_n) / len(i_n)
        I0 = np.fft.fft(i_0) / len(i_0)
        Vp = Vp[:len(positive_freqs)]
        Vn = Vn[:len(positive_freqs)]
        V0 = V0[:len(positive_freqs)]
        Ip = Ip[:len(positive_freqs)]
        In = In[:len(positive_freqs)]
        I0 = I0[:len(positive_freqs)]
        Yp0=Ip[wd]/V0[wd]
        Yn0=In[wd]/V0[wd]
        Y00=I0[wd]/V0[wd]
        Y0pn = np.array([[Ypp, Ypn],
                         [Ynp, Ynn]])
        Y[fd0[i]]=Y0pn
        [Y_mag[fd0[i]], Y_rad[fd0[i]]]=cartz2pol(Y0pn)
    # Open loop for varying the parameter value in the scanner
    # Close PSCAD once the scan has finished
    pscad.quit()    
    # Load the space-separated file into a pandas DataFrame with regex for delimiter
    df3 = pd.read_csv('Frequency_Scan_LAGR.gf81/Ymatrix_pn_unb.csv', delimiter=',', header=None, engine='python')
    # Renombrar las columnas de acuerdo a su significado
    df3.columns = ['jw_str', 'Ym_pp', 'Ya_pp', 'Ym_pn', 'Ya_pn', 'Ym_np', 'Ya_np', 'Ym_nn', 'Ya_nn']
    # Convertir las columnas de strings complejas a números complejos
    jw_sampling = df3['jw_str'].apply(lambda x: x.real).to_numpy()  # Extraer solo la parte real de fw
    Ym_pp = df3['Ym_pp'].to_numpy()
    Ya_pp = df3['Ya_pp'].to_numpy()
    Ym_pn = df3['Ym_pn'].to_numpy()
    Ya_pn = df3['Ya_pn'].to_numpy()
    Ym_np = df3['Ym_np'].to_numpy()
    Ya_np = df3['Ya_np'].to_numpy()
    Ym_nn = df3['Ym_nn'].to_numpy()
    Ya_nn = df3['Ya_nn'].to_numpy()
    fsampling = np.real(jw_sampling)/(2*np.pi)
    # pn0_plot_response(fsampling, Ym_pp, Ya_pp, Ym_pn, Ya_pn, Ym_np, Ya_np, Ym_nn, Ya_nn, fd0, Y_mag, Y_rad)
    pn0_plot(fsampling, fd0, Y_mag, Y_rad)
# Tiempo total de ejecución
end_time = time.time()
execution_time = end_time - start_time
print(f"Tiempo de ejecución: {execution_time} segundos")
    # sys.exit()