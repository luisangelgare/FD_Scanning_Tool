% Technical University of Catalonia (UPC)
% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
% Doctoral Program in Electrical Engineering
% Developed by: Luis Angel Garcia Reyes, MSc

% Frequency domain scanning tool (FDs) for modern power systems 

clear all;
clc;

% Simulation parameters
program='FDScanningTool.slx'; % Program selection
model='FDScanningTool';
Tobs=3; % Observation time
delta_t=10E-6; % Fixed step time
jw1=1i*2*pi*(logspace(0,log10(1/delta_t),10000)); % Complex frequency vector
fd0=unique(round(logspace(-1,log10(1/delta_t)-0.4,10)));  % Perturbation frequencies in Hz
dist_time=1; % Disturbance time
t_w1=1.8; % Starting time window
t_w2=2.8; % Ending time window
Tobs0=t_w2-t_w1; % Total time window
fs=1/Tobs0; % Sampling frequency for FTT
samples_window=round((t_w2-t_w1)/delta_t); % Ensure samples is an integer
time_vector=(t_w1:delta_t:t_w2-delta_t)'; % Time vector for the disturbance
% Fundamental and disturbance sources data
f0=50; % Fundamental base frequency
w0=2*pi*f0; % Angular velocity (rad/s)
Vpeak=1; % Fundamental peak voltage
Ipeak=1; % Fundamental peak current
Rsource=1E-3; % Fundamental series resistance
R_load=1E6; % Termination for transmision lines and others
phi=(2/3)*pi; % Phase of 120º between phases
Vdist_value=0.01*Vpeak; % Voltage disturbance value
Idist_value=0.05*Ipeak; % Current disturbance value

% Scanner settings:
% 1 -> Voltage perturbation 
% 2 -> Current perturbation 
scanner_type=1;

% 1 -> Single-tone perturbation 
% 2 -> RBS perturbation
% 3 -> Multi-tone perturbation
signal_type=1;

% 1 -> ABC scan 
% 2 -> qd0 scan 
% 3 -> pn0 scan
scanner_selector=3;

% 1 -> linear response exists
% 0 -> no linear response is available
linear=1;

% 1 -> RLC load
% 2 -> PI section
% 3 -> Frequency Dependent line
circuit_type=2; 

% 1 -> Balanced 
% 2 -> Unbalanced
balanced_type=2; 

% Signal type building
switch signal_type
    case 1
        % No actions
    case 2
        % RBS or PRBS signals perturbation strategy
        % Voltage perturbation
        V_rbs_ex = idinput([samples_window, 3], 'rbs', [], [0, 0.1]); % RBS generation signals
        Vsignal_dist1 = [time_vector, V_rbs_ex(:, 1)]; % 2-dimension disturbance vector in q [time,rbs signal] 
        Vsignal_dist2 = [time_vector, V_rbs_ex(:, 2)]; % 2-dimension disturbance vector in d [time,rbs signal]
        Vsignal_dist3 = [time_vector, V_rbs_ex(:, 3)]; % 2-dimension disturbance vector in d [time,rbs signal]
        % Current perturbation
        I_rbs_ex = idinput([samples_window, 3], 'prbs', [], [-Idist_value, Idist_value]); % RBS generation signals
        Isignal_dist1 = [time_vector, I_rbs_ex(:, 1)]; % 2-dimension disturbance vector in q [time,rbs signal] 
        Isignal_dist2 = [time_vector, I_rbs_ex(:, 2)]; % 2-dimension disturbance vector in d [time,rbs signal]
        Isignal_dist3 = [time_vector, I_rbs_ex(:, 3)]; % 2-dimension disturbance vector in d [time,rbs signal]
    case 3
        % Cosenoidal signals perturbation strategy
        freq_multiples = [1, 2, 3, 4, 5]; % Multiples of the base frequency
        specific_freqs = fd0; % If non-empty, this overrides freq_multiples
        Vmag = Vdist_value*ones(1,length(specific_freqs)); % Amplitudes for each frequency
        Imag = Idist_value*ones(1,length(specific_freqs)); % Amplitudes for each frequency
        % Combine into a 2-dimensional vector [time, signal]
        Vdisturbance_signal = multisine([fd0(1),fd0(end)], round(1/delta_t), samples_window, ...
              'PhaseResponse', 'Schroeder', ... % Usar fases de Schroeder
              'Normalise', true, ...            % Normalizar la señal
              'StartAtZero', true);             % Iniciar en un cruce por cero
        Idisturbance_signal=Vdisturbance_signal;
        Vsignal_dist1 = [time_vector, 0.1*Vdisturbance_signal']; % a-signal
        Vsignal_dist2=Vsignal_dist1; % b-signal
        Vsignal_dist3=Vsignal_dist1; % c-signal
        Isignal_dist1 = [time_vector, 0.1*Idisturbance_signal']; % a-signal
        Isignal_dist2=Isignal_dist1; % b-signal
        Isignal_dist3=Isignal_dist1; % c-signal
end

% Validation elements available
if linear == 1
% Linear model to validate or black-box parameters (if they exist)
switch circuit_type
    %% RLC load
    case 1
        % RLC circuit data
        R=1; % Resistance
        L=0.1; % Inductance
        C=10E-6; % Capacitance
        switch scanner_selector
            case 1
                % ABC frequency response
            switch balanced_type
                case 1
                    for n=1:length(imag(jw1))
                        Z_RL=R*eye(3)+jw1(n)*L*eye(3);
                        Z_C=inv((jw1(n)*C*eye(3)));
                        Yphases(:,:,n)=inv(Z_RL+Z_C);
                        Zphases(:,:,n)=Z_RL+Z_C;
                    end
                case 2
                    Rdesb=[R    0     0
                               0  0.9*R   0
                               0    0   0.8*R];
                    Ldesb=[L    0     0
                           0  0.9*L   0
                           0    0   0.8*L];
                    Cdesb=[C    0     0
                           0  0.9*C   0
                           0    0   0.8*C];
                    for n=1:length(imag(jw1))
                        Z_RL=Rdesb+jw1(n)*Ldesb;
                        Z_C=inv((jw1(n)*Cdesb));
                        Yphases(:,:,n)=inv(Z_RL+Z_C);
                        Zphases(:,:,n)=Z_RL+Z_C;
                    end
            end
            case 2
                % dq0 frequency response
                Ass=[-R/L  -w0   -1/L    0
                      w0   -R/L   0    -1/L
                      1/C   0     0     -w0
                      0     1/C   w0     0];
                Bss=[1/L  0
                      0  1/L
                      0   0
                      0   0];
                Css=[1 0 0 0
                     0 1 0 0
                     0 0 1 0
                     0 0 0 1];
                Dss=[0 0
                     0 0
                     0 0
                     0 0];
                RLC_state={'X1q','X1d','X2q','Xq2'};
                input_RLC={'Vq','Vd'};
                output_RLC={'Iq','Id','Vcq','Vcd'};
                H_RLC=ss(Ass,Bss,Css,Dss,'StateName',RLC_state,...
                    'inputname',input_RLC,'outputname',output_RLC);
                [Ym_RLC,Ya_RLC]=bode(H_RLC,imag(jw1));
                if scanner_type==2
                [Zm_RLC, Za_RLC]=Zmatrix_qd_ss(length(imag(jw1)), Ym_RLC, Ya_RLC);
                end
            case 3
                % Positive, negative, cero sequences
                alfa=exp(1i*(2/3)*pi);
                Ainv=(1/3)*[1    1     1
                         1  alfa  alfa^2
                         1 alfa^2 alfa];
                A=[1    1      1
                   1  alfa^2  alfa
                   1  alfa    alfa^2];
                switch balanced_type
                    case 1
                % Balanced RLC
                Zrlc_A=R+jw1*L+1./(jw1*C); % RLC frequency response
                Zrlc_B=R+jw1*L+1./(jw1*C); % RLC frequency response
                Zrlc_C=R+jw1*L+1./(jw1*C); % RLC frequency response
                    case 2
                % Unbalanced RLC
                Zrlc_A=R+jw1*L+1./(jw1*C); % RLC frequency response
                Zrlc_B=R*0.9+jw1*L*0.9+1./(jw1*C*0.9); % RLC frequency response
                Zrlc_C=R*0.8+jw1*L*0.8+1./(jw1*C*0.8); % RLC frequency response
                end
                for n=1:length(jw1) 
                    Y_rlc=[1/Zrlc_A(n) 0 0
                           0 1/Zrlc_B(n) 0
                           0 0 1/Zrlc_C(n)];
                    Z_rlc=[Zrlc_A(n) 0 0
                           0 Zrlc_B(n) 0
                           0 0 Zrlc_C(n)];
                    Y0pn=Ainv*Y_rlc*A;
                    Z0pn=Ainv*Z_rlc*A;
                    Yphases(:,:,n)=Y0pn;
                    Zphases(:,:,n)=Z0pn;
                    Ym_RLC(:,:,n)=abs(Y0pn);
                    Ya_RLC(:,:,n)=(180/pi)*angle(Y0pn);
                    Zm_RLC(:,:,n)=abs(Z0pn);
                    Za_RLC(:,:,n)=(180/pi)*angle(Z0pn);
                end
%                 Ymatrix_pn=[transpose(imag(jw1)) ...
%                 squeeze(Ym_RLC(1,1,:)) squeeze(Ya_RLC(1,1,:))...
%                 squeeze(Ym_RLC(1,2,:)) squeeze(Ya_RLC(1,2,:))...
%                 squeeze(Ym_RLC(2,1,:)) squeeze(Ya_RLC(2,1,:))...
%                 squeeze(Ym_RLC(2,2,:)) squeeze(Ya_RLC(2,2,:))];
%                 writematrix(Ymatrix_pn, 'Ymatrix_pn_unb.csv');
        end
    %% PI section
    case 2
        % RLC data
        R=0.1; % Resistance (ohm)
        L=0.001; % Inductance (H)
        C=0.001; % Capacitance (F)
        C1=C/2;
        C2=C1;
        R_load=1E-3;
        switch scanner_selector
            case 1
                % State space representation of PI section
            case 2
                % % State space representation of PI section
                A_source=[0 0
                      0 0];
                B_source=[0 0 0 0
                      0 0 0 0];
                C_source=[0 0
                      0 0];
                D_source=[1/Rsource 0 -1/Rsource 0
                        0    1/Rsource 0 -1/Rsource];
                Xsource_s={'Xq1','Xd2'};
                X_source_in={'vq_source','vd_source' 'v1q','v1d'};
                X_source_out={'iq_source','id_source'};
                H_source=ss(A_source,B_source,C_source,D_source,'StateName',...
                    Xsource_s,'inputname',X_source_in,'outputname',X_source_out);
                % C1 shunt state space
                A_C1=[0      -w0
                      w0       0];
                B_C1=[1/C1     0      -1/C1      0
                        0       1/C1       0      -1/C1];
                C_C1=[1 0
                      0 1];
                D_C1=[0 0 0 0
                      0 0 0 0];
                X_C1_s={'Vq_C1','Vd_C1'};
                X_C1_in={'iq_source','id_source','iq_RL','id_RL'};
                X_C1_out={'v1q','v1d'};
                H_C1=ss(A_C1,B_C1,C_C1,D_C1,'StateName',...
                    X_C1_s,'inputname',X_C1_in,'outputname',X_C1_out);
                % RL branch state space
                A_RLb=[-R/L       -w0
                          w0        -R/L];
                B_RLb=[1/L     0    -1/L       0
                        0     1/L     0       -1/L];
                C_RLb=[1 0
                      0 1];
                D_RLb=[0 0 0 0
                      0 0 0 0]; 
                X_RL_s={'Iq_RL','Id_RL'};
                X_RL_in={'v1q','v1d','v2q' 'v2d'};
                X_RL_out={'iq_RL','id_RL'};
                H_RLb=ss(A_RLb,B_RLb,C_RLb,D_RLb,'StateName',...
                    X_RL_s,'inputname',X_RL_in,'outputname',X_RL_out);
                % C2 shunt state space
                A_C2=[0      -w0
                      w0       0];
                B_C2=[1/C2     0      -1/C2      0
                        0       1/C2       0      -1/C2];
                C_C2=[1 0
                      0 1];
                D_C2=[0 0 0 0
                      0 0 0 0];
                X_C2_s={'Vq_C2','Vd_C2'};
                X_C2_in={'iq_RL','id_RL','iq_load','id_load'};
                X_C2_out={'v2q' 'v2d'};
                H_C2=ss(A_C2,B_C2,C_C2,D_C2,'StateName',...
                    X_C2_s,'inputname',X_C2_in,'outputname',X_C2_out);
                % Rload state space
                A_R=[0 0
                     0 0];
                B_R=[0 0
                     0 0];
                C_R=[0 0
                     0 0];
                D_R=[1/R_load      0    
                         0      1/R_load];
                X_R_s={'Iq_R','Id_R'};
                X_R_in={'v2q','v2d'};
                X_R_out={'iq_load','id_load'};
                H_R=ss(A_R,B_R,C_R,D_R,'StateName',...
                    X_R_s,'inputname',X_R_in,'outputname',X_R_out);
                % Trasmission line admitance
                inPI={'vq_source','vd_source'};
                outPI={'iq_source','id_source'};
                H_PI_TL=connect(H_source,H_C1,H_RLb,H_C2,H_R,inPI,outPI);
                % RLC frequency response
                % jw1=2*pi*logspace(0,4,2^14);
                [Ym_RLC,Ya_RLC]=bode(H_PI_TL,imag(jw1));
                if scanner_type==2
                [Zm_RLC, Za_RLC]=Zmatrix_qd_ss(length(imag(jw1)), Ym_RLC, Ya_RLC);
                end
            case 3
                R_load=1E-3;
                % Positive, negative, cero sequences
                alfa=exp(1i*(2/3)*pi);
                Ainv=(1/3)*[1    1     1
                1  alfa  alfa^2
                1 alfa^2 alfa];
                A=[1    1      1
                1  alfa^2  alfa
                1  alfa    alfa^2];
                switch balanced_type
                    case 1
                % Balanced RLC
                ZA=R+jw1*L; % Series impedance
                YA=jw1*(C/2); % Shunt admittance
                ZB=R+jw1*L; % Series impedance
                YB=jw1*(C/2); % Shunt admittance
                ZC=R+jw1*L; % Series impedance
                YC=jw1*(C/2); % Shunt admittance
                case 2
                % Unbalanced RLC
                ZA=R+jw1*L; % Series impedance
                YA=jw1*(C/2); % Shunt admittance
                ZB=R*0.9+jw1*L*0.9; % Series impedance
                YB=jw1*(C/2)*0.9; % Shunt admittance
                ZC=R*0.8+jw1*L*0.8; % Series impedance
                YC=jw1*(C/2)*0.8; % Shunt admittance
                end
                % % Frequency domain impedance and admitance   
                for n=1:length(imag(jw1))
                Yabc1=[YA(n)+1/ZA(n)       0           0
                          0          YB(n)+1/ZB(n)     0
                          0                0        YC(n)+1/ZC(n)];
                Yabc0(:,:,n)=Yabc1;
                Y0pn=Ainv*Yabc1*A;
                Y_pn(:,:,n)=Y0pn(2:3,2:3);
                Ym_RLC(:,:,n)=abs(Y0pn);
                Ya_RLC(:,:,n)=(180/pi)*angle(Y0pn);
                Z0pn=inv(Y0pn);
                Z_pn(:,:,n)=Z0pn(2:3,2:3);
                Zm_RLC(:,:,n)=abs(Z0pn);
                Za_RLC(:,:,n)=(180/pi)*angle(Z0pn);
                end
        end
    %% FD cable or overhead line
    case 3
        switch scanner_selector
            case 1
                R_load=1E-6;
                % Positive, negative, cero sequences
                alfa=exp(1i*(2/3)*pi);
                Ainv=(1/3)*[1    1     1
                            1  alfa  alfa^2
                            1 alfa^2 alfa];
                A=[1    1      1
                   1  alfa^2  alfa
                   1  alfa    alfa^2];
                s=1E-8+jw1; % Vector de frecuencia compleja s
                % Ro=28.3E-9; % Resistividad del aluminio en Ohm/m
                % Ro=5.8341e-8;
                Ro=3.206E-8;
                Rsuelo=100; % Resistividad del terreno en Ohm/m
                SigmaSuelo=1/Rsuelo; % Conductividad del terreno en S/m
                % SigmaAl=37E6; % Conductividad del aluminio en S/m
                % SigmaAl=1/4.021E-8;
                SigmaAl=1/Ro;
                SigmaAire=8E-15; % Conductividad del aire en S/m
                % SigmaAire=0; % Conductividad del aire en S/m
                Epsilon=8.854187E-12; % Permitividad en el vacío en F/m
                Miu=4*pi*1E-7; % Permeabilidad en el vacío en Wb/A.m
                % Linea 
                numconducfase=1; % Número de conductores por fase
                numfases=3; % Número de fases
                numguardas=2; % Número de conductores de guarda
                lineas=5; % Número de lineas del sistema (3fases+2guardas=5lineas)
                tr=3.5; % Cooeficiente de trenzado (1.5 ≤ t ≤ 3.5)
                R=(3.09E-2)/2; % Radio del haz de conductores
                Mxy=[-5.49       14.4   R    % x1,y1,r1
                      0          15.62  R    % x2,y2,r2
                      5.49       14.4   R    % x3,y3,r3
                     -3.05       18.21  (1.26E-2)/2   % xG1,yG1,rG1
                      3.05       18.21  (1.26E-2)/2]; % xG2,yG2,rG2
%                 [Ap,Bp,Z_abc,Y_abc]=ParametrosAB(s,numfases,lineas,numconducfase,Miu,...
%                     SigmaSuelo,SigmaAl,SigmaAire,Epsilon,Ro,tr,100E3,R,Mxy);
                [Z_abc_TL,Y_abc_TL]=LineParameters(s);
                [Ap,Bp]=Nodal_param(Z_abc_TL,Y_abc_TL,100E3);
                % Admitance matrix in frequency domain
                for n=1:length(imag(jw1))
                Yp=(squeeze(Ap(:,:,n))-squeeze(Bp(:,:,n)));
                Ys=(squeeze(Bp(:,:,n)));
%                 Yabc1=inv(inv(Yp+Ys)+inv(Yp+ones(3)*(1/R_load)));
                Yabc1=(Yp+Ys);
                % Y_pn(:,:,n)=ReduccionKron(Y0pn,1);
                Yphase(:,n)=(Yabc1(1,1));
                Yphases(:,:,n)=Yabc1;
                Zphases(:,:,n)=inv(Yabc1);
                end
            case 2
                R_load=1E-3;
            case 3
                R_load=1E-6;
                % Positive, negative, cero sequences
                alfa=exp(1i*(2/3)*pi);
                Ainv=(1/3)*[1    1     1
                            1  alfa  alfa^2
                            1 alfa^2 alfa];
                A=[1    1      1
                   1  alfa^2  alfa
                   1  alfa    alfa^2];
%                 s=1E-8+jw1; % Vector de frecuencia compleja s
%                 Ro=28.3E-9; % Resistividad del aluminio en Ohm/m
%                 Ro=5.8341e-8;
                Ro=3.206E-8;
                Rsuelo=100; % Resistividad del terreno en Ohm/m
                SigmaSuelo=1/Rsuelo; % Conductividad del terreno en S/m
                % SigmaAl=37E6; % Conductividad del aluminio en S/m
                % SigmaAl=1/4.021E-8;
                SigmaAl=1/Ro;
                SigmaAire=0; % Conductividad del aire en S/m
                % SigmaAire=0; % Conductividad del aire en S/m
                Epsilon=8.854187E-12; % Permitividad en el vacío en F/m
                Miu=4*pi*1E-7; % Permeabilidad en el vacío en Wb/A.m
                % Linea 
                numconducfase=1; % Número de conductores por fase
                numfases=3; % Número de fases
                numguardas=2; % Número de conductores de guarda
                lineas=5; % Número de lineas del sistema (3fases+2guardas=5lineas)
                tr=3.5; % Cooeficiente de trenzado (1.5 ≤ t ≤ 3.5)
                R=(3.556E-2)/2; % Radio del haz de conductores
                % Matriz de coordenadas en el plano XY con radios de los conductores
                Mxy=[-12.8016      20.7265   R    % x1,y1,r1
                      0            20.7265   R    % x2,y2,r2
                      12.8016      20.7265   R    % x3,y3,r3
                     -8.9916       32.9185  (1.27E-2)/2   % xG1,yG1,rG1
                      8.9916       32.9185  (1.27E-2)/2]; % xG2,yG2,rG2
%                 [Ap,Bp,Z_abc,Y_abc]=ParametrosAB(jw1,numfases,lineas,numconducfase,Miu,...
%                     SigmaSuelo,SigmaAl,SigmaAire,Epsilon,Ro,tr,100E3,R,Mxy);
                [Z_abc_TL,Y_abc_TL]=LineParameters(jw1);
                [Ap,Bp]=Nodal_param(Z_abc_TL,Y_abc_TL,100E3);
                % Admitance matrix in frequency domain  
                for n=1:length(imag(jw1))
                Yp=(squeeze(Ap(:,:,n))-squeeze(Bp(:,:,n)));
                Ys=(squeeze(Bp(:,:,n)));
                Yabc1=(Yp+Ys);
                Yabc0(:,:,n)=Yabc1;
                Y0pn=Ainv*Yabc1*A;
                Z0pn=inv(Y0pn);
                Y0pn_full(:,:,n)=Y0pn;
                Z0pn_full(:,:,n)=Z0pn;
                Y_pn=Y0pn(2:3,2:3);
                Z_pn=Z0pn(2:3,2:3);
                Ym_RLC(:,:,n)=abs(Y_pn);
                Ya_RLC(:,:,n)=(180/pi)*angle(Y_pn);
                Zm_RLC(:,:,n)=abs(Z_pn);
                Za_RLC(:,:,n)=(180/pi)*angle(Z_pn);
%                 Ymatrix_pn=[transpose(imag(jw1)) ...
%                 squeeze(Ym_RLC(1,1,:)) squeeze(Ya_RLC(1,1,:))...
%                 squeeze(Ym_RLC(1,2,:)) squeeze(Ya_RLC(1,2,:))...
%                 squeeze(Ym_RLC(2,1,:)) squeeze(Ya_RLC(2,1,:))...
%                 squeeze(Ym_RLC(2,2,:)) squeeze(Ya_RLC(2,2,:))];
%                 writematrix(Ymatrix_pn, 'Ymatrix_pn_unb.csv');
                end
        end
end
end
 
%% Comprehensive frequency scan tool methodology

tic % clock start

% Subsystem enabling
voltage_type=[model, '/Frequency Scanner/Voltage Strategy']; 
current_type=[model, '/Frequency Scanner/Current Strategy'];
ABC_V_scanner=[model, '/Frequency Scanner/Voltage Strategy/ABCdisturbance']; 
qd0_V_scanner=[model, '/Frequency Scanner/Voltage Strategy/qd0disturbance'];
pn0_V_scanner=[model, '/Frequency Scanner/Voltage Strategy/pn0disturbance'];
ABC_I_scanner=[model, '/Frequency Scanner/Current Strategy/ABCdisturbance2']; 
qd0_I_scanner=[model, '/Frequency Scanner/Current Strategy/qd0disturbance2'];
pn0_I_scanner=[model, '/Frequency Scanner/Current Strategy/pn0disturbance2'];
multi_tone_V_ABC=[model, '/Frequency Scanner/Voltage Strategy/ABCdisturbance/Multi-tone'];
single_tone_V_ABC=[model, '/Frequency Scanner/Voltage Strategy/ABCdisturbance/Single-tone'];
multi_tone_I_ABC=[model, '/Frequency Scanner/Current Strategy/ABCdisturbance2/Multi-tone'];
single_tone_I_ABC=[model, '/Frequency Scanner/Current Strategy/ABCdisturbance2/Single-tone'];

switch scanner_selector
    %% ABC Scanner
    case 1
        set_param(ABC_V_scanner, 'Commented', 'off');
        set_param(qd0_V_scanner, 'Commented', 'on');
        set_param(pn0_V_scanner, 'Commented', 'on');
        set_param(ABC_I_scanner, 'Commented', 'off');
        set_param(qd0_I_scanner, 'Commented', 'on');
        set_param(pn0_I_scanner, 'Commented', 'on');
        if scanner_type==1
            set_param(voltage_type, 'Commented', 'off');
            set_param(current_type, 'Commented', 'on');
            if signal_type==1
                set_param(multi_tone_V_ABC, 'Commented', 'on');
                set_param(single_tone_V_ABC, 'Commented', 'off');
                for n=1:length(fd0)
                    fd=fd0(n);
                    % a-injection
                    dist_value_a=Vdist_value;
                    dist_value_b=0;
                    dist_value_c=0;
                    out=sim(program);
                    td1=find((out.tout)>=t_w1,1); % t1 time window
                    td2=find((out.tout)>=t_w2,1); % t2 time window
                    % Vector windowed extraction
                    va=out.Vabc(td1:td2,1);
                    vb=out.Vabc(td1:td2,2);
                    vc=out.Vabc(td1:td2,3);
                    ia=out.Iabc(td1:td2,1);
                    ib=out.Iabc(td1:td2,2);
                    ic=out.Iabc(td1:td2,3);
                    Va=fft(va)/length(va);
                    Vb=fft(vb)/length(vb);
                    Vc=fft(vc)/length(vc);
                    Ia=fft(ia)/length(ia);
                    Ib=fft(ib)/length(ib);
                    Ic=fft(ic)/length(ic);
                    wd=round(fd/fs)+1;
                    Yaa=Ia(wd)/Va(wd);
                    Yba=Ib(wd)/Va(wd);
                    Yca=Ic(wd)/Va(wd);
                    % b-injection
                    dist_value_a=0;
                    dist_value_b=Vdist_value;
                    dist_value_c=0;
                    out=sim(program);
                    va=out.Vabc(td1:td2,1);
                    vb=out.Vabc(td1:td2,2);
                    vc=out.Vabc(td1:td2,3);
                    ia=out.Iabc(td1:td2,1);
                    ib=out.Iabc(td1:td2,2);
                    ic=out.Iabc(td1:td2,3);
                    Va=fft(va)/length(va);
                    Vb=fft(vb)/length(vb);
                    Vc=fft(vc)/length(vc);
                    Ia=fft(ia)/length(ia);
                    Ib=fft(ib)/length(ib);
                    Ic=fft(ic)/length(ic);
                    Yab=Ia(wd)/Vb(wd);
                    Ybb=Ib(wd)/Vb(wd);
                    Ycb=Ic(wd)/Vb(wd);
                    % c-injection
                    dist_value_a=0;
                    dist_value_b=0;
                    dist_value_c=Vdist_value;
                    out=sim(program);
                    va=out.Vabc(td1:td2,1);
                    vb=out.Vabc(td1:td2,2);
                    vc=out.Vabc(td1:td2,3);
                    ia=out.Iabc(td1:td2,1);
                    ib=out.Iabc(td1:td2,2);
                    ic=out.Iabc(td1:td2,3);
                    Va=fft(va)/length(va);
                    Vb=fft(vb)/length(vb);
                    Vc=fft(vc)/length(vc);
                    Ia=fft(ia)/length(ia);
                    Ib=fft(ib)/length(ib);
                    Ic=fft(ic)/length(ic);
                    Yac=Ia(wd)/Vc(wd);
                    Ybc=Ib(wd)/Vc(wd);
                    Ycc=Ic(wd)/Vc(wd);
                    Yabc=[Yaa Yab Yac
                          Yba Ybb Ybc
                          Yca Ycb Ycc];
                    Y_abc(:,:,n)=Yabc;
                end
                Ya=squeeze(Y_abc(1,1,:));
                Yb=squeeze(Y_abc(2,2,:));
                Yc=squeeze(Y_abc(3,3,:));
            else
                set_param(multi_tone_V_ABC, 'Commented', 'off');
                set_param(single_tone_V_ABC, 'Commented', 'on');
                % a-injection
                dist_value_a=Vsignal_dist1;
                dist_value_b=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                dist_value_c=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                % Vector windowed extraction
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va=fft(va)/length(va);
                Vb=fft(vb)/length(vb);
                Vc=fft(vc)/length(vc);
                Ia=fft(ia)/length(ia);
                Ib=fft(ib)/length(ib);
                Ic=fft(ic)/length(ic);
                Yaa=Ia./Va;
                Yba=Ib./Va;
                Yca=Ic./Va;
                % b-injection
                dist_value_b=Vsignal_dist2;
                dist_value_a=[Vsignal_dist2(:,1),zeros(samples_window,1)];
                dist_value_c=[Vsignal_dist2(:,1),zeros(samples_window,1)];
                out=sim(program);
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va=fft(va)/length(va);
                Vb=fft(vb)/length(vb);
                Vc=fft(vc)/length(vc);
                Ia=fft(ia)/length(ia);
                Ib=fft(ib)/length(ib);
                Ic=fft(ic)/length(ic);
                Yab=Ia./Vb;
                Ybb=Ib./Vb;
                Ycb=Ic./Vb;
                % c-injection
                dist_value_c=Vsignal_dist3;
                dist_value_a=[Vsignal_dist3(:,1),zeros(samples_window,1)];
                dist_value_b=[Vsignal_dist3(:,1),zeros(samples_window,1)];
                out=sim(program);
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va=fft(va)/length(va);
                Vb=fft(vb)/length(vb);
                Vc=fft(vc)/length(vc);
                Ia=fft(ia)/length(ia);
                Ib=fft(ib)/length(ib);
                Ic=fft(ic)/length(ic);
                Yac=Ia./Vc;
                Ybc=Ib./Vc;
                Ycc=Ic./Vc;
                for n=1:length(Ycc)
                    Y_abc(:,:,n)=[Yaa(n) Yab(n) Yac(n)
                                  Yba(n) Ybb(n) Ybc(n)
                                  Yca(n) Ycb(n) Ycc(n)];
                end
                fd0=0:fs:1/delta_t;
            end
        else
            set_param(voltage_type, 'Commented', 'on');
            set_param(current_type, 'Commented', 'off');
            if signal_type==1
                set_param(multi_tone_I_ABC, 'Commented', 'on');
                set_param(single_tone_I_ABC, 'Commented', 'off');
            for n=1:length(fd0)
                fd=fd0(n);
                % a-injection
                dist_value_a=Idist_value;
                dist_value_b=0;
                dist_value_c=0;
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                % Vector windowed extraction
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va_a=fft(va)/length(va);
                Vb_a=fft(vb)/length(vb);
                Vc_a=fft(vc)/length(vc);
                Ia_a=fft(ia)/length(ia);
                Ib_a=fft(ib)/length(ib);
                Ic_a=fft(ic)/length(ic);
                % b-injection
                dist_value_a=0;
                dist_value_b=Idist_value;
                dist_value_c=0;
                out=sim(program);
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va_b=fft(va)/length(va);
                Vb_b=fft(vb)/length(vb);
                Vc_b=fft(vc)/length(vc);
                Ia_b=fft(ia)/length(ia);
                Ib_b=fft(ib)/length(ib);
                Ic_b=fft(ic)/length(ic);
                % c-injection
                dist_value_a=0;
                dist_value_b=0;
                dist_value_c=Idist_value;
                out=sim(program);
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va_c=fft(va)/length(va);
                Vb_c=fft(vb)/length(vb);
                Vc_c=fft(vc)/length(vc);
                Ia_c=fft(ia)/length(ia);
                Ib_c=fft(ib)/length(ib);
                Ic_c=fft(ic)/length(ic);
                wd=round(fd/fs)+1;
                Vabc=[Va_a(wd) Va_b(wd) Va_c(wd)
                      Vb_a(wd) Vb_b(wd) Vb_c(wd)
                      Vc_a(wd) Vc_b(wd) Vc_c(wd)];
                Iabc=[Ia_a(wd) Ia_b(wd) Ia_c(wd)
                      Ib_a(wd) Ib_b(wd) Ib_c(wd)
                      Ic_a(wd) Ic_b(wd) Ic_c(wd)];
                Z_abc(:,:,n)=Vabc*inv(Iabc);
            end
            Za=squeeze(Z_abc(1,1,:));
            Zb=squeeze(Z_abc(2,2,:));
            Zc=squeeze(Z_abc(3,3,:));
            else
                set_param(multi_tone_I_ABC, 'Commented', 'off');
                set_param(single_tone_I_ABC, 'Commented', 'on');
                % a-injection
                dist_value_a=Isignal_dist1;
                dist_value_b=[Isignal_dist1(:,1),zeros(samples_window,1)];
                dist_value_c=[Isignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                % Vector windowed extraction
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va_a=fft(va)/length(va);
                Vb_a=fft(vb)/length(vb);
                Vc_a=fft(vc)/length(vc);
                Ia_a=fft(ia)/length(ia);
                Ib_a=fft(ib)/length(ib);
                Ic_a=fft(ic)/length(ic);
                % b-injection
                dist_value_b=Isignal_dist2;
                dist_value_a=[Isignal_dist2(:,1),zeros(samples_window,1)];
                dist_value_c=[Isignal_dist2(:,1),zeros(samples_window,1)];
                out=sim(program);
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va_b=fft(va)/length(va);
                Vb_b=fft(vb)/length(vb);
                Vc_b=fft(vc)/length(vc);
                Ia_b=fft(ia)/length(ia);
                Ib_b=fft(ib)/length(ib);
                Ic_b=fft(ic)/length(ic);
                % c-injection
                dist_value_c=Isignal_dist3;
                dist_value_a=[Isignal_dist3(:,1),zeros(samples_window,1)];
                dist_value_b=[Isignal_dist3(:,1),zeros(samples_window,1)];
                out=sim(program);
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va_c=fft(va)/length(va);
                Vb_c=fft(vb)/length(vb);
                Vc_c=fft(vc)/length(vc);
                Ia_c=fft(ia)/length(ia);
                Ib_c=fft(ib)/length(ib);
                Ic_c=fft(ic)/length(ic);
            for n=1:length(Ic_c)
                Vabc=[Va_a(n) Va_b(n) Va_c(n)
                      Vb_a(n) Vb_b(n) Vb_c(n)
                      Vc_a(n) Vc_b(n) Vc_c(n)];
                Iabc=[Ia_a(n) Ia_b(n) Ia_c(n)
                      Ib_a(n) Ib_b(n) Ib_c(n)
                      Ic_a(n) Ic_b(n) Ic_c(n)];
                Z_abc(:,:,n)=Vabc*inv(Iabc);
            end
            fd0=0:fs:1/delta_t;
            Za=squeeze(Z_abc(1,1,:));
            Zb=squeeze(Z_abc(2,2,:));
            Zc=squeeze(Z_abc(3,3,:));
            end
        end
        switch linear
            case 1
                if scanner_type==1
                    ABCPlot(fd0, Y_abc, Yphases, jw1);
                else
                    ABCPlot(fd0, Z_abc, Zphases, jw1);
                end
            case 0
                if scanner_type==1
                    ABCPlot2(fd0, Y_abc);
                else
                    ABCPlot2(fd0, Z_abc);
                end
        end
    %% qd0 Scanner
    case 2
        set_param(ABC_V_scanner, 'Commented', 'on');
        set_param(qd0_V_scanner, 'Commented', 'off');
        set_param(pn0_V_scanner, 'Commented', 'on');
        set_param(ABC_I_scanner, 'Commented', 'on');
        set_param(qd0_I_scanner, 'Commented', 'off');
        set_param(pn0_I_scanner, 'Commented', 'on');
        if scanner_type==1
            set_param(voltage_type, 'Commented', 'off');
            set_param(current_type, 'Commented', 'on');
%             if signal_type==1
                    for n=1:length(fd0)
                        fd=fd0(n); % Frequency value
                        clear dist_value_q dist_value_d act
                        dist_value_q=Vdist_value; % q disturbance
                        dist_value_d=0.0; % d disturbance
                        out=sim(program);
                        tss=find((out.tout)>=dist_time-0.1,1); % t1 time window
                        td1=find((out.tout)>=t_w1,1); % t1 time window
                        td2=find((out.tout)>=t_w2,1); % t2 time window
                        % Vector windowed extraction
                        vq_ss=out.Vqd(tss,1);
                        vd_ss=out.Vqd(tss,2);
                        iq_ss=out.Iqd(tss,1);
                        id_ss=out.Iqd(tss,2);
                        vq=out.Vqd(td1:td2,1);
                        vd=out.Vqd(td1:td2,2);
                        iq=out.Iqd(td1:td2,1);
                        id=out.Iqd(td1:td2,2);
                        % FFT of the time-windowed signals
                        Vq=fft(vq-vq_ss)/length(vq);
                        Iq=fft(iq-iq_ss)/length(iq);
                        Id=fft(id-id_ss)/length(id);
                        wd=round(fd/fs)+1;
                        % dq and qq Admitance calculation in FD
                        Ydq(:,n)=Id(wd)/Vq(wd);
                        Yqq(:,n)=Iq(wd)/Vq(wd);
                        clear dist_value_q dist_value_d act
                        dist_value_q=0.0; % q disturbance
                        dist_value_d=Vdist_value; % d disturbance
                        out=sim(program);
                        td1=find((out.tout)>=t_w1,1);
                        td2=find((out.tout)>=t_w2,1);
                        % Vector windowed extraction
                        vq_ss=out.Vqd(tss,1);
                        vd_ss=out.Vqd(tss,2);
                        iq_ss=out.Iqd(tss,1);
                        id_ss=out.Iqd(tss,2);
                        vq=out.Vqd(td1:td2,1);
                        vd=out.Vqd(td1:td2,2);
                        iq=out.Iqd(td1:td2,1);
                        id=out.Iqd(td1:td2,2);
                        % FFT of the time-windowed signals
                        Vd=fft(vd-vd_ss)/length(vd);
                        Iq=fft(iq-iq_ss)/length(iq);
                        Id=fft(id-id_ss)/length(id);
                         % qd and dd Admitance calculation in FD
                        Yqd(:,n)=Iq(wd)/Vd(wd);
                        Ydd(:,n)=Id(wd)/Vd(wd);
                    end
%             else
%                 clear Vdist_signal_q Vdist_signal_d
%                 Vdist_signal_q=Vsignal_dist1;
%                 Vdist_signal_d=[Vsignal_dist1(:,1),zeros(samples_window,1)];
%                 out=sim(program);
%                 tss=find((out.tout)>=dist_time-0.1,1); % t1 time window
%                 td1=find((out.tout)>=t_w1,1); % t1 time window
%                 td2=find((out.tout)>=t_w2,1); % t2 time window
%                 % Vector windowed extraction
%                 vq_ss=out.Vqd(tss,1);
%                 vd_ss=out.Vqd(tss,2);
%                 iq_ss=out.Iqd(tss,1);
%                 id_ss=out.Iqd(tss,2);
%                 vq=out.Vqd(td1:td2,1);
%                 vd=out.Vqd(td1:td2,2);
%                 iq=out.Iqd(td1:td2,1);
%                 id=out.Iqd(td1:td2,2);
%                 % FFT of the time-windowed signals
%                 Vq=fft(vq-vq_ss)/length(vq);
%                 Iq=fft(iq-iq_ss)/length(iq);
%                 Id=fft(id-id_ss)/length(id);
%                 % dq and qq Admitance calculation in FD
%                 Ydq=Id./Vq;
%                 Yqq=Iq./Vq;
%                 clear Vdist_signal_q Vdist_signal_d
%                 Vdist_signal_q=[Vsignal_dist2(:,1),zeros(samples_window,1)];
%                 Vdist_signal_d=Vsignal_dist2;
%                 out=sim(program);
%                 td1=find((out.tout)>=t_w1,1);
%                 td2=find((out.tout)>=t_w2,1);
%                 % Vector windowed extraction
%                 vq_ss=out.Vqd(tss,1);
%                 vd_ss=out.Vqd(tss,2);
%                 iq_ss=out.Iqd(tss,1);
%                 id_ss=out.Iqd(tss,2);
%                 vq=out.Vqd(td1:td2,1);
%                 vd=out.Vqd(td1:td2,2);
%                 iq=out.Iqd(td1:td2,1);
%                 id=out.Iqd(td1:td2,2);
%                 % FFT of the time-windowed signals
%                 Vd=fft(vd-vd_ss)/length(vd);
%                 Iq=fft(iq-iq_ss)/length(iq);
%                 Id=fft(id-id_ss)/length(id);
%                  % qd and dd Admitance calculation in FD
%                 Yqd=Iq./Vd;
%                 Ydd=Id./Vd;
%             end
        else
            set_param(voltage_type, 'Commented', 'on');
            set_param(current_type, 'Commented', 'off');
            for n=1:length(fd0)
                fd=fd0(n); % Frequency value
                clear dist_value_q dist_value_d act
                dist_value_q=Idist_value; % q disturbance
                dist_value_d=0.0; % d disturbance
                out=sim(program);
                tss=find((out.tout)>=dist_time-0.1,1); % t1 time window
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                % Vector windowed extraction
                vq_ss=out.Vqd(tss,1);
                vd_ss=out.Vqd(tss,2);
                iq_ss=out.Iqd(tss,1);
                id_ss=out.Iqd(tss,2);
                vq=out.Vqd(td1:td2,1);
                vd=out.Vqd(td1:td2,2);
                iq=out.Iqd(td1:td2,1);
                id=out.Iqd(td1:td2,2);
                % FFT of the time-windowed signals
                Vq_q=fft(vq-vq_ss)/length(vq);
                Vd_q=fft(vd-vd_ss)/length(vd);
                Iq_q=fft(iq-iq_ss)/length(iq);
                Id_q=fft(id-id_ss)/length(id);
                % d-injection
                clear dist_value_q dist_value_d act
                dist_value_q=0.0; % q disturbance
                dist_value_d=Idist_value; % d disturbance
                out=sim(program);
                td1=find((out.tout)>=t_w1,1);
                td2=find((out.tout)>=t_w2,1);
                % Vector windowed extraction
                vq_ss=out.Vqd(tss,1);
                vd_ss=out.Vqd(tss,2);
                iq_ss=out.Iqd(tss,1);
                id_ss=out.Iqd(tss,2);
                vq=out.Vqd(td1:td2,1);
                vd=out.Vqd(td1:td2,2);
                iq=out.Iqd(td1:td2,1);
                id=out.Iqd(td1:td2,2);
                % FFT of the time-windowed signals
                Vq_d=fft(vq-vq_ss)/length(vq);
                Vd_d=fft(vd-vd_ss)/length(vd);
                Iq_d=fft(iq-iq_ss)/length(iq);
                Id_d=fft(id-id_ss)/length(id);
                 % qd and dd Admitance calculation in FD
                 wd=round(fd/fs)+1;
                 Vqd=[Vq_q(wd) Vq_d(wd)
                      Vd_q(wd) Vd_d(wd)];
                 Iqd=[Iq_q(wd) Iq_d(wd)
                      Id_q(wd) Id_d(wd)];
                Z_qd=Vqd*inv(Iqd);
                Zqq(:,n)=Z_qd(1,1);
                Zqd(:,n)=Z_qd(1,2);
                Zdq(:,n)=Z_qd(2,1);
                Zdd(:,n)=Z_qd(2,2);
            end
        end
            switch linear
                case 1
                    if scanner_type==1
                    qd0Plot(fd0, jw1, Ym_RLC, Ya_RLC, Yqq, Yqd, Ydq, Ydd);
                    else
                    qd0Plot(fd0, jw1, Zm_RLC, Za_RLC, Zqq, Zqd, Zdq, Zdd);
                    end
                case 0
                    if scanner_type==1
                    qd0Plot2(fd0, Yqq, Yqd, Ydq, Ydd);
                    else
                    qd0Plot2(fd0, Zqq, Zqd, Zdq, Zdd);
                    end
            end
    %% pn0 Scanner
    case 3
        set_param(ABC_V_scanner, 'Commented', 'on');
        set_param(qd0_V_scanner, 'Commented', 'on');
        set_param(pn0_V_scanner, 'Commented', 'off');
        set_param(ABC_I_scanner, 'Commented', 'on');
        set_param(qd0_I_scanner, 'Commented', 'on');
        set_param(pn0_I_scanner, 'Commented', 'off');
        if scanner_type==1
            set_param(voltage_type, 'Commented', 'off');
            set_param(current_type, 'Commented', 'on');
            for n=1:length(fd0)
                fd=fd0(n); % Frequency value
                clear dist_value_p dist_value_n dist_value_0
                dist_value_0=0.0; % d disturbance
                dist_value_p=Vdist_value; % q disturbance
                dist_value_n=0.0; % d disturbance
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                time0=out.tout(td1:td2,1); %
                % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0=fft(v0)/length(v0);
                Vp=fft(vp)/length(vp);
                Vn=fft(vn)/length(vn);
                I0=fft(i0)/length(i0);
                Ip=fft(ip)/length(ip);
                In=fft(in)/length(in);
                wd=round(fd/fs)+1;
                % dq and qq Admitance calculation in FD
                Y_0p0=I0(wd)/Vp(wd);
                Y_pp0=Ip(wd)/Vp(wd);
                Y_np0=In(wd)/Vp(wd);
                Ypp(:,n)=Y_pp0;
                Ynp(:,n)=Y_np0;
                clear dist_value_p dist_value_n dist_value_0
                dist_value_0=0.0; % d disturbance
                dist_value_p=0.0; % q disturbance
                dist_value_n=Vdist_value; % d disturbance
                fd=-fd0(n);
                out=sim(program);
                % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0=fft(v0)/length(v0);
                Vp=fft(vp)/length(vp);
                Vn=fft(vn)/length(vn);
                I0=fft(i0)/length(i0);
                Ip=fft(ip)/length(ip);
                In=fft(in)/length(in);
                Y_0n0=I0(wd)/Vn(wd);
                Y_pn0=Ip(wd)/Vn(wd);
                Y_nn0=In(wd)/Vn(wd);
                Ypn(:,n)=Y_pn0;
                Ynn(:,n)=Y_nn0;
                clear dist_value_p dist_value_n dist_value_0
                dist_value_0=Vdist_value; % d disturbance
                dist_value_p=0.0; % q disturbance
                dist_value_n=0.0; % d disturbance
                fd=fd0(n);
                out=sim(program);
                % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0=fft(v0)/length(v0);
                Vp=fft(vp)/length(vp);
                Vn=fft(vn)/length(vn);
                I0=fft(i0)/length(i0);
                Ip=fft(ip)/length(ip);
                In=fft(in)/length(in);
                 % qd and dd Admitance calculation in FD
                Y_000=I0(wd)/V0(wd);
                Y_p00=Ip(wd)/V0(wd);
                Y_n00=In(wd)/V0(wd);
                Y_0pn=[Y_000 Y_0p0 Y_0n0
                       Y_p00 Y_pp0 Y_pn0
                       Y_n00 Y_np0 Y_nn0];
                Y0pn_all(:,:,n)=Y_0pn;
            end
        else
            set_param(voltage_type, 'Commented', 'on');
            set_param(current_type, 'Commented', 'off');
            for n=1:length(fd0)
                fd=fd0(n); % Frequency value
                clear dist_value_p dist_value_n dist_value_0
                dist_value_0=0.0; % d disturbance
                dist_value_p=Idist_value; % q disturbance
                dist_value_n=0.0; % d disturbance
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                time0=out.tout(td1:td2,1); %
                % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0_p=fft(v0)/length(v0);
                Vp_p=fft(vp)/length(vp);
                Vn_p=fft(vn)/length(vn);
                I0_p=fft(i0)/length(i0);
                Ip_p=fft(ip)/length(ip);
                In_p=fft(in)/length(in);
                clear dist_value_p dist_value_n dist_value_0
                dist_value_0=0.0; % d disturbance
                dist_value_p=0.0; % q disturbance
                dist_value_n=Idist_value; % d disturbance
                fd=-fd0(n);
                out=sim(program);
                % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0_n=fft(v0)/length(v0);
                Vp_n=fft(vp)/length(vp);
                Vn_n=fft(vn)/length(vn);
                I0_n=fft(i0)/length(i0);
                Ip_n=fft(ip)/length(ip);
                In_n=fft(in)/length(in);
                clear dist_value_p dist_value_n dist_value_0
                dist_value_0=Idist_value; % d disturbance
                dist_value_p=0.0; % q disturbance
                dist_value_n=0.0; % d disturbance
                fd=fd0(n);
                out=sim(program);
                % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0_0=fft(v0)/length(v0);
                Vp_0=fft(vp)/length(vp);
                Vn_0=fft(vn)/length(vn);
                I0_0=fft(i0)/length(i0);
                Ip_0=fft(ip)/length(ip);
                In_0=fft(in)/length(in);
                 % 0pn and dd Admitance calculation in FD
                wd=round(fd/fs)+1;
                V0pn=[V0_0(wd) V0_p(wd) V0_n(wd)
                      Vp_0(wd) Vp_p(wd) Vp_n(wd)
                      Vn_0(wd) Vn_p(wd) Vn_n(wd)];
                I0pn=[I0_0(wd) I0_p(wd) I0_n(wd)
                      Ip_0(wd) Ip_p(wd) Ip_n(wd)
                      In_0(wd) In_p(wd) In_n(wd)];
                Z_0pn=V0pn*inv(I0pn);
                Zpp(:,n)=Z_0pn(2,2);
                Zpn(:,n)=Z_0pn(2,3);
                Znp(:,n)=Z_0pn(3,2);
                Znn(:,n)=Z_0pn(3,3);
                Z0pn_all(:,:,n)=Z_0pn;
            end
         end
        switch linear
            case 1
                if scanner_type==1
                    pn0Plot(fd0, jw1, Ym_RLC, Ya_RLC, Ypp, Ypn, Ynp, Ynn);
                else
                    pn0Plot(fd0, jw1, Zm_RLC, Za_RLC, Zpp, Zpn, Znp, Znn); 
                end
            case 0
                if scanner_type==1
                    pn0Plot2(fd0, Ypp, Ypn, Ynp, Ynn);
                else
                    pn0Plot2(fd0, Zpp, Zpn, Znp, Znn); 
                end
        end
end
toc