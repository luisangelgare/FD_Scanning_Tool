%% Technical University of Catalonia (UPC)
%% Higher Technical School of Industrial Engineering of Barcelona (ETSEIB)
%% Centre of Technological Innovation in Static Converters and Drives (CITCEA)
%% Doctoral Program in Electrical Engineering
%% Developed by: Luis Angel Garcia Reyes, MSc
%% Frequency domain scanning tool (FDs) for modern power systems 

% GNU General Public License v3.0 (GPL-3.0)
% Copyright (C) 2025 Luis Angel Garcia Reyes, UPC-MSCA-ADOreD
% Email: luis.reyes@upc.edu
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see <https://www.gnu.org/licenses/>.

% This work has received funding from the ADOreD project
% under the European Union’s Horizon Europe Research and 
% Innovation Programme under the Marie Skłodowska-Curie 
% Grant Agreement No. 101073554.

%% Preliminar parameters based on the simulation settings

fsampling=round(1/delta_t); % Sampling frequency for the simulation 
Tobs0=1/fs; % Total time window
t_w1=Tinit; % Starting time window
t_w2=Tinit+Tobs0; % Ending time window
Tobs=t_w2; % Full observation time
dist_time_ss=Tobs+10000; % No disturbance for calculating the steady state
dist_time_f=Tinit-0.3; % effective disturbance time for the scanner
time_vector=(t_w1:delta_t:t_w2)'; % Time vector for the disturbance
samples_window=length(time_vector); % Ensure samples is an integer
Rsource=1E-6; % Fundamental series resistance (suggested to avoid errors)
phi=(2/3)*pi; % Phase of 120º between phases
Vdist_value=0.01*Vpeak; % Voltage disturbance value
Idist_value=0.03*Ipeak; % Current disturbance value
dist_time=dist_time_ss; % Set the first disturbance time for steady state
% Logic to determine the ss_type
if Vq_ss == 0 && Vd_ss == 0
    Vss_type = 1; % Three-phase ss voltage source
else
    Vss_type = 2; % dq0 to ABC ss voltage source
end
if Iq_ss == 0 && Id_ss == 0
    Iss_type = 1; % Three-phase ss current source
else
    Iss_type = 2; % dq0 to ABC ss current source
end

%% Signal disturbance development with multi-sine and random binary strategies

disp('SIaD Tool is starting...')
disp('---')
fprintf('Simulation set with %s perturbation, %s signal, %s frame.\n', ...
    choose(scanner_type, 'voltage', 'current'), ...
    choose(signal_type, 'single-tone', 'PRBS', 'multi-tone'), ...
    choose(scanner_selector, 'ABC', 'dq0', '0pn'));

switch signal_type
    case 1
        % No actions
    case 2
        % RBS or PRBS signals perturbation strategy
        % Voltage perturbation
        V_rbs_ex = idinput([samples_window, 3], 'prbs', [], [0, Vdist_value]); % RBS generation signals
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
        freq_multiples = [1, 2, 3, 4, 5, 6, 7, 8]; % Multiples of the base frequency
        specific_freqs = fd0; % If non-empty, this overrides freq_multiples
        Vmag = Vdist_value*ones(1,length(specific_freqs)); % Amplitudes for each frequency
        Imag = Idist_value*ones(1,length(specific_freqs)); % Amplitudes for each frequency
        % Combine into a 2-dimensional vector [time, signal]
        Vdisturbance_signal = multisine([fd0(1),fd0(end)], fsampling, samples_window, ...
              'PhaseResponse', 'Schroeder', ... %  Schroeder phases
              'Normalise', true, ...            % Norm of the signal
              'StartAtZero', true);             % Initialize in a zero-crossing point
        Idisturbance_signal=Vdisturbance_signal;
        Vsignal_dist1 = [time_vector, Vdist_value*(Vdisturbance_signal)']; % a-signal
        Vsignal_dist2=Vsignal_dist1; % b-signal
        Vsignal_dist3=Vsignal_dist1; % c-signal
        Isignal_dist1 = [time_vector, Idist_value*Idisturbance_signal']; % a-signal
        Isignal_dist2=Isignal_dist1; % b-signal
        Isignal_dist3=Isignal_dist1; % c-signal
end

%% Subsystem subblocks identification and assignment

tic % clock start
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
multi_tone_V_qd0=[model, '/Frequency Scanner/Voltage Strategy/qd0disturbance/Multi-tone'];
single_tone_V_qd0=[model, '/Frequency Scanner/Voltage Strategy/qd0disturbance/Single-tone'];
multi_tone_I_qd0=[model, '/Frequency Scanner/Current Strategy/qd0disturbance2/Multi-tone'];
single_tone_I_qd0=[model, '/Frequency Scanner/Current Strategy/qd0disturbance2/Single-tone'];
multi_tone_V_0pn=[model, '/Frequency Scanner/Voltage Strategy/pn0disturbance/Multi-tone'];
single_tone_V_0pn=[model, '/Frequency Scanner/Voltage Strategy/pn0disturbance/Single-tone'];
multi_tone_I_0pn=[model, '/Frequency Scanner/Current Strategy/pn0disturbance2/Multi-tone'];
single_tone_I_0pn=[model, '/Frequency Scanner/Current Strategy/pn0disturbance2/Single-tone'];

%% Comprehensive frequency scanning tool methodology (abs, dq0, pn0)

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
            fd=0; % Frequency value
            if signal_type==1
                set_param(multi_tone_V_ABC, 'Commented', 'on');
                set_param(single_tone_V_ABC, 'Commented', 'off');
                dist_value_a=0;
                dist_value_b=0;
                dist_value_c=0;
            else
                set_param(multi_tone_V_ABC, 'Commented', 'off');
                set_param(single_tone_V_ABC, 'Commented', 'on');
                dist_value_a=zeros(samples_window,2);
                dist_value_b=zeros(samples_window,2);
                dist_value_c=zeros(samples_window,2);
            end
                disp('SIaD Tool is obtaining the steady state...')
                disp('---')
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                va_ss=out.Vabc(td1:td2,1);
                vb_ss=out.Vabc(td1:td2,2);
                vc_ss=out.Vabc(td1:td2,3);
                ia_ss=out.Iabc(td1:td2,1);
                ib_ss=out.Iabc(td1:td2,2);
                ic_ss=out.Iabc(td1:td2,3);
                dist_time=dist_time_f;
                disp('SIaD Tool is running the system identification process...')
                disp('---')
            if signal_type==1
                for n=1:length(fd0)
                    fd=fd0(n);
                    % a-injection
                    dist_value_a=Vdist_value;
                    dist_value_b=0;
                    dist_value_c=0;
                    out=sim(program);
                    % Vector windowed extraction
                    va=out.Vabc(td1:td2,1);
                    vb=out.Vabc(td1:td2,2);
                    vc=out.Vabc(td1:td2,3);
                    ia=out.Iabc(td1:td2,1);
                    ib=out.Iabc(td1:td2,2);
                    ic=out.Iabc(td1:td2,3);
                    Va=fft(va-va_ss)/length(va);
                    Vb=fft(vb-vb_ss)/length(vb);
                    Vc=fft(vc-vc_ss)/length(vc);
                    Ia=fft(ia-ia_ss)/length(ia);
                    Ib=fft(ib-ib_ss)/length(ib);
                    Ic=fft(ic-ic_ss)/length(ic);
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
                    Va=fft(va-va_ss)/length(va);
                    Vb=fft(vb-vb_ss)/length(vb);
                    Vc=fft(vc-vc_ss)/length(vc);
                    Ia=fft(ia-ia_ss)/length(ia);
                    Ib=fft(ib-ib_ss)/length(ib);
                    Ic=fft(ic-ic_ss)/length(ic);
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
                    Va=fft(va-va_ss)/length(va);
                    Vb=fft(vb-vb_ss)/length(vb);
                    Vc=fft(vc-vc_ss)/length(vc);
                    Ia=fft(ia-ia_ss)/length(ia);
                    Ib=fft(ib-ib_ss)/length(ib);
                    Ic=fft(ic-ic_ss)/length(ic);
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
                % a-injection
                dist_value_a=Vsignal_dist1;
                dist_value_b=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                dist_value_c=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                % Vector windowed extraction
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va=fft(va-va_ss)/length(va);
                Vb=fft(vb-vb_ss)/length(vb);
                Vc=fft(vc-vc_ss)/length(vc);
                Ia=fft(ia-ia_ss)/length(ia);
                Ib=fft(ib-ib_ss)/length(ib);
                Ic=fft(ic-ic_ss)/length(ic);
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
                Va=fft(va-va_ss)/length(va);
                Vb=fft(vb-vb_ss)/length(vb);
                Vc=fft(vc-vc_ss)/length(vc);
                Ia=fft(ia-ia_ss)/length(ia);
                Ib=fft(ib-ib_ss)/length(ib);
                Ic=fft(ic-ic_ss)/length(ic);
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
                Va=fft(va-va_ss)/length(va);
                Vb=fft(vb-vb_ss)/length(vb);
                Vc=fft(vc-vc_ss)/length(vc);
                Ia=fft(ia-ia_ss)/length(ia);
                Ib=fft(ib-ib_ss)/length(ib);
                Ic=fft(ic-ic_ss)/length(ic);
                Yac=Ia./Vc;
                Ybc=Ib./Vc;
                Ycc=Ic./Vc;
                for n=1:length(Ycc)
                    Y_abc(:,:,n)=[Yaa(n) Yab(n) Yac(n)
                                  Yba(n) Ybb(n) Ybc(n)
                                  Yca(n) Ycb(n) Ycc(n)];
                end
                fd0=(1:samples_window)*fs;
            end
        else
            set_param(voltage_type, 'Commented', 'on');
            set_param(current_type, 'Commented', 'off');
            fd=0; % Frequency value
            if signal_type==1
                set_param(multi_tone_I_ABC, 'Commented', 'on');
                set_param(single_tone_I_ABC, 'Commented', 'off');
                dist_value_a=0;
                dist_value_b=0;
                dist_value_c=0;
            else
                set_param(multi_tone_I_ABC, 'Commented', 'off');
                set_param(single_tone_I_ABC, 'Commented', 'on');
                dist_value_a=zeros(samples_window,2);
                dist_value_b=zeros(samples_window,2);
                dist_value_c=zeros(samples_window,2);
            end
                disp('SIaD Tool is obtaining the steady state...')
                disp('---')    
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                va_ss=out.Vabc(td1:td2,1);
                vb_ss=out.Vabc(td1:td2,2);
                vc_ss=out.Vabc(td1:td2,3);
                ia_ss=out.Iabc(td1:td2,1);
                ib_ss=out.Iabc(td1:td2,2);
                ic_ss=out.Iabc(td1:td2,3);
                dist_time=dist_time_f;
                disp('SIaD Tool is running the system identification process...')
                disp('---')
            if signal_type==1
                for n=1:length(fd0)
                    fd=fd0(n);
                    % a-injection
                    dist_value_a=Idist_value;
                    dist_value_b=0;
                    dist_value_c=0;
                    out=sim(program);
                    % Vector windowed extraction
                    va=out.Vabc(td1:td2,1);
                    vb=out.Vabc(td1:td2,2);
                    vc=out.Vabc(td1:td2,3);
                    ia=out.Iabc(td1:td2,1);
                    ib=out.Iabc(td1:td2,2);
                    ic=out.Iabc(td1:td2,3);
                    Va_a=fft(va-va_ss)/length(va);
                    Vb_a=fft(vb-vb_ss)/length(vb);
                    Vc_a=fft(vc-vc_ss)/length(vc);
                    Ia_a=fft(ia-ia_ss)/length(ia);
                    Ib_a=fft(ib-ib_ss)/length(ib);
                    Ic_a=fft(ic-ic_ss)/length(ic);
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
                    Va_b=fft(va-va_ss)/length(va);
                    Vb_b=fft(vb-vb_ss)/length(vb);
                    Vc_b=fft(vc-vc_ss)/length(vc);
                    Ia_b=fft(ia-ia_ss)/length(ia);
                    Ib_b=fft(ib-ib_ss)/length(ib);
                    Ic_b=fft(ic-ic_ss)/length(ic);
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
                    Va_c=fft(va-va_ss)/length(va);
                    Vb_c=fft(vb-vb_ss)/length(vb);
                    Vc_c=fft(vc-vc_ss)/length(vc);
                    Ia_c=fft(ia-ia_ss)/length(ia);
                    Ib_c=fft(ib-ib_ss)/length(ib);
                    Ic_c=fft(ic-ic_ss)/length(ic);
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
                % a-injection
                dist_value_a=Isignal_dist1;
                dist_value_b=[Isignal_dist1(:,1),zeros(samples_window,1)];
                dist_value_c=[Isignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                % Vector windowed extraction
                va=out.Vabc(td1:td2,1);
                vb=out.Vabc(td1:td2,2);
                vc=out.Vabc(td1:td2,3);
                ia=out.Iabc(td1:td2,1);
                ib=out.Iabc(td1:td2,2);
                ic=out.Iabc(td1:td2,3);
                Va_a=fft(va-va_ss)/length(va);
                Vb_a=fft(vb-vb_ss)/length(vb);
                Vc_a=fft(vc-vc_ss)/length(vc);
                Ia_a=fft(ia-ia_ss)/length(ia);
                Ib_a=fft(ib-ib_ss)/length(ib);
                Ic_a=fft(ic-ic_ss)/length(ic);
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
                Va_b=fft(va-va_ss)/length(va);
                Vb_b=fft(vb-vb_ss)/length(vb);
                Vc_b=fft(vc-vc_ss)/length(vc);
                Ia_b=fft(ia-ia_ss)/length(ia);
                Ib_b=fft(ib-ib_ss)/length(ib);
                Ic_b=fft(ic-ic_ss)/length(ic);
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
                Va_c=fft(va-va_ss)/length(va);
                Vb_c=fft(vb-vb_ss)/length(vb);
                Vc_c=fft(vc-vc_ss)/length(vc);
                Ia_c=fft(ia-ia_ss)/length(ia);
                Ib_c=fft(ib-ib_ss)/length(ib);
                Ic_c=fft(ic-ic_ss)/length(ic);
            for n=1:length(Ic_c)
                Vabc=[Va_a(n) Va_b(n) Va_c(n)
                      Vb_a(n) Vb_b(n) Vb_c(n)
                      Vc_a(n) Vc_b(n) Vc_c(n)];
                Iabc=[Ia_a(n) Ia_b(n) Ia_c(n)
                      Ib_a(n) Ib_b(n) Ib_c(n)
                      Ic_a(n) Ic_b(n) Ic_c(n)];
                Z_abc(:,:,n)=Vabc*inv(Iabc);
            end
            fd0=(1:samples_window)*fs;
            Za=squeeze(Z_abc(1,1,:));
            Zb=squeeze(Z_abc(2,2,:));
            Zc=squeeze(Z_abc(3,3,:));
            end
        end
        switch linear
            case 1
                if scanner_type==1
                    ABCPlot(fd0, Y_abc, Yphases, jw1);
                    fprintf('SIaD Tool finished. Results stored in Ya, Yb and Yc.\n');
                else
                    ABCPlot(fd0, Z_abc, Zphases, jw1);
                    fprintf('SIaD Tool finished. Results stored in Za, Zb and Zc.\n');
                end
            case 0
                if scanner_type==1
                    ABCPlot2(fd0, Y_abc);
                    fprintf('SIaD Tool finished. Results stored in Ya, Yb and Yc.\n');
                else
                    ABCPlot2(fd0, Z_abc);
                    fprintf('SIaD Tool finished. Results stored in Za, Zb and Zc.\n');
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
            fd=0; % Frequency value
                if signal_type==1
                    set_param(multi_tone_V_qd0, 'Commented', 'on');
                    set_param(single_tone_V_qd0, 'Commented', 'off');
                    dist_value_q=0.0; % q disturbance
                    dist_value_d=0.0; % d disturbance
                else
                    set_param(multi_tone_V_qd0, 'Commented', 'off');
                    set_param(single_tone_V_qd0, 'Commented', 'on');
                    dist_value_q=zeros(samples_window,2);
                    dist_value_d=zeros(samples_window,2);
                end
                disp('SIaD Tool is obtaining the steady state...')
                disp('---')
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                vq_ss=out.Vqd(td1:td2,1);
                vd_ss=out.Vqd(td1:td2,2);
                iq_ss=out.Iqd(td1:td2,1);
                id_ss=out.Iqd(td1:td2,2);
                dist_time=dist_time_f;
                disp('SIaD Tool is running the system identification process...')
                disp('---')
                if signal_type==1
                            for n=1:length(fd0)
                                fd=fd0(n); % Frequency value
                                clear dist_value_q dist_value_d act
                                dist_value_q=Vdist_value; % q disturbance
                                dist_value_d=0.0; % d disturbance
                                out=sim(program);
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
                else
                    % q-injection
                    clear dist_value_q dist_value_d act
                    dist_value_q=Vsignal_dist1;
                    dist_value_d=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                    out=sim(program);
                    % Vector windowed extraction
                    vq=out.Vqd(td1:td2,1);
                    vd=out.Vqd(td1:td2,2);
                    iq=out.Iqd(td1:td2,1);
                    id=out.Iqd(td1:td2,2);
                    Vq=fft(vq-vq_ss)/length(vq);
                    Vd=fft(vd-vd_ss)/length(vd);
                    Iq=fft(iq-iq_ss)/length(iq);
                    Id=fft(id-id_ss)/length(id);
                    Yqq=Iq./Vq;
                    Ydq=Id./Vq;
                    % d-injection
                    clear dist_value_q dist_value_d act
                    dist_value_d=Vsignal_dist1;
                    dist_value_q=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                    out=sim(program);
                    % Vector windowed extraction
                    vq=out.Vqd(td1:td2,1);
                    vd=out.Vqd(td1:td2,2);
                    iq=out.Iqd(td1:td2,1);
                    id=out.Iqd(td1:td2,2);
                    Vq=fft(vq-vq_ss)/length(vq);
                    Vd=fft(vd-vd_ss)/length(vd);
                    Iq=fft(iq-iq_ss)/length(iq);
                    Id=fft(id-id_ss)/length(id);
                    Ydd=Id./Vd;
                    Yqd=Iq./Vd;
                    % for n=1:length(Yqq)
                    %     Y_abc(:,:,n)=[Yaa(n) Yab(n) Yac(n)
                    %                   Yba(n) Ybb(n) Ybc(n)
                    %                   Yca(n) Ycb(n) Ycc(n)];
                    % end
                    fd0=(1:samples_window)*fs;
                end
        else
            set_param(voltage_type, 'Commented', 'on');
            set_param(current_type, 'Commented', 'off');
            fd=0; % Frequency value
            if signal_type==1
                set_param(multi_tone_I_qd0, 'Commented', 'on');
                set_param(single_tone_I_qd0, 'Commented', 'off');
                dist_value_q=0.0; % q disturbance
                dist_value_d=0.0; % d disturbance
            else
                set_param(multi_tone_I_qd0, 'Commented', 'off');
                set_param(single_tone_I_qd0, 'Commented', 'on');
                dist_value_q=zeros(samples_window,2);
                dist_value_d=zeros(samples_window,2);
            end
            disp('SIaD Tool is obtaining the steady state...')
            disp('---')
            out=sim(program);
            td1=find((out.tout)>=t_w1,1); % t1 time window
            td2=find((out.tout)>=t_w2,1); % t2 time window
            vq_ss=out.Vqd(td1:td2,1);
            vd_ss=out.Vqd(td1:td2,2);
            iq_ss=out.Iqd(td1:td2,1);
            id_ss=out.Iqd(td1:td2,2);
            dist_time=dist_time_f;
            disp('SIaD Tool is running the system identification process...')
                disp('---')
            if signal_type==1
                for n=1:length(fd0)
                    fd=fd0(n);
                    % q-injection
                    dist_value_q=Idist_value;
                    dist_value_d=0;
                    out=sim(program);
                    % Vector windowed extraction
                    vq=out.Vqd(td1:td2,1);
                    vd=out.Vqd(td1:td2,2);
                    iq=out.Iqd(td1:td2,1);
                    id=out.Iqd(td1:td2,2);
                    Vq_q=fft(vq-vq_ss)/length(vq);
                    Vd_q=fft(vd-vd_ss)/length(vd);
                    Iq_q=fft(iq-iq_ss)/length(iq);
                    Id_q=fft(id-id_ss)/length(id);
                    % d-injection
                    dist_value_d=Idist_value;
                    dist_value_q=0;
                    out=sim(program);
                    % Vector windowed extraction
                    vq=out.Vqd(td1:td2,1);
                    vd=out.Vqd(td1:td2,2);
                    iq=out.Iqd(td1:td2,1);
                    id=out.Iqd(td1:td2,2);
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
            else
                % q-injection
                dist_value_q=Isignal_dist1;
                dist_value_d=[Isignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                % Vector windowed extraction
                vq=out.Vqd(td1:td2,1);
                vd=out.Vqd(td1:td2,2);
                iq=out.Iqd(td1:td2,1);
                id=out.Iqd(td1:td2,2);
                Vq_q=fft(vq-vq_ss)/length(vq);
                Vd_q=fft(vd-vd_ss)/length(vd);
                Iq_q=fft(iq-iq_ss)/length(iq);
                Id_q=fft(id-id_ss)/length(id);
                % d-injection
                dist_value_d=Isignal_dist1;
                dist_value_q=[Isignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                % Vector windowed extraction
                vq=out.Vqd(td1:td2,1);
                vd=out.Vqd(td1:td2,2);
                iq=out.Iqd(td1:td2,1);
                id=out.Iqd(td1:td2,2);
                Vq_d=fft(vq-vq_ss)/length(vq);
                Vd_d=fft(vd-vd_ss)/length(vd);
                Iq_d=fft(iq-iq_ss)/length(iq);
                Id_d=fft(id-id_ss)/length(id);
            for n=1:length(Id_d)
                Vqd=[Vq_q(n) Vq_d(n)
                     Vd_q(n) Vd_d(n)];
                Iqd=[Iq_q(n) Iq_d(n)
                     Id_q(n) Id_d(n)];
                Z_qd=Vqd*inv(Iqd);
                Zqq(:,n)=Z_qd(1,1);
                Zqd(:,n)=Z_qd(1,2);
                Zdq(:,n)=Z_qd(2,1);
                Zdd(:,n)=Z_qd(2,2);
            end
            fd0=(1:samples_window)*fs;
            end
        end
            switch linear
                case 1
                    if scanner_type==1
                    qd0Plot(fd0, jw1, Ym_RLC, Ya_RLC, Yqq, Yqd, Ydq, Ydd);
                    fprintf('SIaD Tool finished. Results stored in Yqq, Yqd, Ydq and Ydd.\n');
                    else
                    qd0Plot(fd0, jw1, Zm_RLC, Za_RLC, Zqq, Zqd, Zdq, Zdd);
                    fprintf('SIaD Tool finished. Results stored in Zqq, Zqd, Zdq and Zdd.\n');
                    end
                case 0
                    if scanner_type==1
                    qd0Plot2(fd0, Yqq, Yqd, Ydq, Ydd);
                    fprintf('SIaD Tool finished. Results stored in Yqq, Yqd, Ydq and Ydd.\n');
                    else
                    qd0Plot2(fd0, Zqq, Zqd, Zdq, Zdd);
                    fprintf('SIaD Tool finished. Results stored in Zqq, Zqd, Zdq and Zdd.\n');
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
            fd=0; % Frequency value
            if signal_type==1
                set_param(multi_tone_V_0pn, 'Commented', 'on');
                set_param(single_tone_V_0pn, 'Commented', 'off');
                dist_value_0=0;
                dist_value_p=0;
                dist_value_n=0;
            else
                set_param(multi_tone_V_0pn, 'Commented', 'off');
                set_param(single_tone_V_0pn, 'Commented', 'on');
                dist_value_0=zeros(samples_window,2);
                dist_value_p=zeros(samples_window,2);
                dist_value_n=zeros(samples_window,2);
            end
                disp('SIaD Tool is obtaining the steady state...')
                disp('---')    
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                v0_ss=out.V0pn(td1:td2,1);
                vp_ss=out.V0pn(td1:td2,2);
                vn_ss=out.V0pn(td1:td2,3);
                i0_ss=out.I0pn(td1:td2,1);
                ip_ss=out.I0pn(td1:td2,2);
                in_ss=out.I0pn(td1:td2,3);
                dist_time=dist_time_f;
                disp('SIaD Tool is running the system identification process...')
                disp('---')
            if signal_type==1
                for n=1:length(fd0)
                    fd=fd0(n);
                    % p-injection
                    dist_value_p=Vdist_value;
                    dist_value_n=0;
                    dist_value_0=0;
                    out=sim(program);
                    % Vector windowed extraction
                    v0=out.V0pn(td1:td2,1);
                    vp=out.V0pn(td1:td2,2);
                    vn=out.V0pn(td1:td2,3);
                    i0=out.I0pn(td1:td2,1);
                    ip=out.I0pn(td1:td2,2);
                    in=out.I0pn(td1:td2,3);
                    V0=fft(v0-v0_ss)/length(v0);
                    Vp=fft(vp-vp_ss)/length(vp);
                    Vn=fft(vn-vn_ss)/length(vn);
                    I0=fft(i0-i0_ss)/length(i0);
                    Ip=fft(ip-ip_ss)/length(ip);
                    In=fft(in-in_ss)/length(in);
                    wd=round(fd/fs)+1;
                    Y_0p0=I0(wd)/Vp(wd);
                    Y_pp0=Ip(wd)/Vp(wd);
                    Y_np0=In(wd)/Vp(wd);
                    Ypp(:,n)=Y_pp0;
                    Ynp(:,n)=Y_np0;
                    % n-injection
                    dist_value_p=0;
                    dist_value_n=Vdist_value;
                    dist_value_0=0;
                    out=sim(program);
                    % Vector windowed extraction
                    v0=out.V0pn(td1:td2,1);
                    vp=out.V0pn(td1:td2,2);
                    vn=out.V0pn(td1:td2,3);
                    i0=out.I0pn(td1:td2,1);
                    ip=out.I0pn(td1:td2,2);
                    in=out.I0pn(td1:td2,3);
                    % FFT of the time-windowed signals
                    V0=fft(v0-v0_ss)/length(v0);
                    Vp=fft(vp-vp_ss)/length(vp);
                    Vn=fft(vn-vn_ss)/length(vn);
                    I0=fft(i0-i0_ss)/length(i0);
                    Ip=fft(ip-ip_ss)/length(ip);
                    In=fft(in-in_ss)/length(in);
                    Y_0n0=I0(wd)/Vn(wd);
                    Y_pn0=Ip(wd)/Vn(wd);
                    Y_nn0=In(wd)/Vn(wd);
                    Ypn(:,n)=Y_pn0;
                    Ynn(:,n)=Y_nn0;
                    % 0-injection
                    dist_value_p=0;
                    dist_value_n=0;
                    dist_value_0=Vdist_value;
                    out=sim(program);
                    % Vector windowed extraction
                    v0=out.V0pn(td1:td2,1);
                    vp=out.V0pn(td1:td2,2);
                    vn=out.V0pn(td1:td2,3);
                    i0=out.I0pn(td1:td2,1);
                    ip=out.I0pn(td1:td2,2);
                    in=out.I0pn(td1:td2,3);
                    % FFT of the time-windowed signals
                    V0=fft(v0-v0_ss)/length(v0);
                    Vp=fft(vp-vp_ss)/length(vp);
                    Vn=fft(vn-vn_ss)/length(vn);
                    I0=fft(i0-i0_ss)/length(i0);
                    Ip=fft(ip-ip_ss)/length(ip);
                    In=fft(in-in_ss)/length(in);
                     % 0pn and dd Admitance calculation in FD
                    Y_000=I0(wd)/V0(wd);
                    Y_p00=Ip(wd)/V0(wd);
                    Y_n00=In(wd)/V0(wd);
                    Y_0pn=[Y_000 Y_0p0 Y_0n0
                           Y_p00 Y_pp0 Y_pn0
                           Y_n00 Y_np0 Y_nn0];
                    Y0pn_all(:,:,n)=Y_0pn;
                end
            else
                % p-injection
                dist_value_p=Vsignal_dist1;
                dist_value_n=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                dist_value_0=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                 % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0=fft(v0-v0_ss)/length(v0);
                Vp=fft(vp-vp_ss)/length(vp);
                Vn=fft(vn-vn_ss)/length(vn);
                I0=fft(i0-i0_ss)/length(i0);
                Ip=fft(ip-ip_ss)/length(ip);
                In=fft(in-in_ss)/length(in);
                Y0p=I0./Vp;
                Ypp=Ip./Vp;
                Ynp=In./Vp;
                % n-injection
                dist_value_n=Vsignal_dist2;
                dist_value_p=[Vsignal_dist2(:,1),zeros(samples_window,1)];
                dist_value_0=[Vsignal_dist2(:,1),zeros(samples_window,1)];
                out=sim(program);
                 % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0=fft(v0-v0_ss)/length(v0);
                Vp=fft(vp-vp_ss)/length(vp);
                Vn=fft(vn-vn_ss)/length(vn);
                I0=fft(i0-i0_ss)/length(i0);
                Ip=fft(ip-ip_ss)/length(ip);
                In=fft(in-in_ss)/length(in);
                Y0n=I0./Vn;
                Ypn=Ip./Vn;
                Ynn=In./Vn;
                % 0-injection
                dist_value_0=Vsignal_dist3;
                dist_value_p=[Vsignal_dist3(:,1),zeros(samples_window,1)];
                dist_value_n=[Vsignal_dist3(:,1),zeros(samples_window,1)];
                out=sim(program);
                 % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0=fft(v0-v0_ss)/length(v0);
                Vp=fft(vp-vp_ss)/length(vp);
                Vn=fft(vn-vn_ss)/length(vn);
                I0=fft(i0-i0_ss)/length(i0);
                Ip=fft(ip-ip_ss)/length(ip);
                In=fft(in-in_ss)/length(in);
                Y00=I0./V0;
                Yp0=Ip./V0;
                Yn0=In./V0;
                for n=1:length(Yn0)
                    Y_0pn(:,:,n)=[Y00(n) Y0p(n) Y0n(n)
                                  Yp0(n) Ypp(n) Ypn(n)
                                  Yn0(n) Ynp(n) Ynn(n)];
                end
                fd0=(1:samples_window)*fs;
            end
        else
            set_param(voltage_type, 'Commented', 'on');
            set_param(current_type, 'Commented', 'off');
            fd=0; % Frequency value
            if signal_type==1
                set_param(multi_tone_I_0pn, 'Commented', 'on');
                set_param(single_tone_I_0pn, 'Commented', 'off');
                dist_value_0=0;
                dist_value_p=0;
                dist_value_n=0;
            else
                set_param(multi_tone_I_0pn, 'Commented', 'off');
                set_param(single_tone_I_0pn, 'Commented', 'on');
                dist_value_0=zeros(samples_window,2);
                dist_value_p=zeros(samples_window,2);
                dist_value_n=zeros(samples_window,2);
            end
                disp('SIaD Tool is obtaining the steady state...')
                disp('---')    
                out=sim(program);
                td1=find((out.tout)>=t_w1,1); % t1 time window
                td2=find((out.tout)>=t_w2,1); % t2 time window
                v0_ss=out.V0pn(td1:td2,1);
                vp_ss=out.V0pn(td1:td2,2);
                vn_ss=out.V0pn(td1:td2,3);
                i0_ss=out.I0pn(td1:td2,1);
                ip_ss=out.I0pn(td1:td2,2);
                in_ss=out.I0pn(td1:td2,3);
                dist_time=dist_time_f;
                disp('SIaD Tool is running the system identification process...')
                disp('---')
            if signal_type==1
                for n=1:length(fd0)
                    fd=fd0(n);
                    % p-injection
                    dist_value_p=Vdist_value;
                    dist_value_n=0;
                    dist_value_0=0;
                    out=sim(program);
                    % Vector windowed extraction
                    v0=out.V0pn(td1:td2,1);
                    vp=out.V0pn(td1:td2,2);
                    vn=out.V0pn(td1:td2,3);
                    i0=out.I0pn(td1:td2,1);
                    ip=out.I0pn(td1:td2,2);
                    in=out.I0pn(td1:td2,3);
                    V0_p=fft(v0-v0_ss)/length(v0);
                    Vp_p=fft(vp-vp_ss)/length(vp);
                    Vn_p=fft(vn-vn_ss)/length(vn);
                    I0_p=fft(i0-i0_ss)/length(i0);
                    Ip_p=fft(ip-ip_ss)/length(ip);
                    In_p=fft(in-in_ss)/length(in);
                    % n-injection
                    dist_value_p=0;
                    dist_value_n=Vdist_value;
                    dist_value_0=0;
                    out=sim(program);
                    % Vector windowed extraction
                    v0=out.V0pn(td1:td2,1);
                    vp=out.V0pn(td1:td2,2);
                    vn=out.V0pn(td1:td2,3);
                    i0=out.I0pn(td1:td2,1);
                    ip=out.I0pn(td1:td2,2);
                    in=out.I0pn(td1:td2,3);
                    % FFT of the time-windowed signals
                    V0_n=fft(v0-v0_ss)/length(v0);
                    Vp_n=fft(vp-vp_ss)/length(vp);
                    Vn_n=fft(vn-vn_ss)/length(vn);
                    I0_n=fft(i0-i0_ss)/length(i0);
                    Ip_n=fft(ip-ip_ss)/length(ip);
                    In_n=fft(in-in_ss)/length(in);
                    % 0-injection
                    dist_value_p=0;
                    dist_value_n=0;
                    dist_value_0=Vdist_value;
                    out=sim(program);
                    % Vector windowed extraction
                    v0=out.V0pn(td1:td2,1);
                    vp=out.V0pn(td1:td2,2);
                    vn=out.V0pn(td1:td2,3);
                    i0=out.I0pn(td1:td2,1);
                    ip=out.I0pn(td1:td2,2);
                    in=out.I0pn(td1:td2,3);
                    % FFT of the time-windowed signals
                    V0_0=fft(v0-v0_ss)/length(v0);
                    Vp_0=fft(vp-vp_ss)/length(vp);
                    Vn_0=fft(vn-vn_ss)/length(vn);
                    I0_0=fft(i0-i0_ss)/length(i0);
                    Ip_0=fft(ip-ip_ss)/length(ip);
                    In_0=fft(in-in_ss)/length(in);
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
            else
                % p-injection
                dist_value_p=Vsignal_dist1;
                dist_value_n=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                dist_value_0=[Vsignal_dist1(:,1),zeros(samples_window,1)];
                out=sim(program);
                 % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0_p=fft(v0-v0_ss)/length(v0);
                Vp_p=fft(vp-vp_ss)/length(vp);
                Vn_p=fft(vn-vn_ss)/length(vn);
                I0_p=fft(i0-i0_ss)/length(i0);
                Ip_p=fft(ip-ip_ss)/length(ip);
                In_p=fft(in-in_ss)/length(in);
                % n-injection
                dist_value_n=Vsignal_dist2;
                dist_value_p=[Vsignal_dist2(:,1),zeros(samples_window,1)];
                dist_value_0=[Vsignal_dist2(:,1),zeros(samples_window,1)];
                out=sim(program);
                 % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0_n=fft(v0-v0_ss)/length(v0);
                Vp_n=fft(vp-vp_ss)/length(vp);
                Vn_n=fft(vn-vn_ss)/length(vn);
                I0_n=fft(i0-i0_ss)/length(i0);
                Ip_n=fft(ip-ip_ss)/length(ip);
                In_n=fft(in-in_ss)/length(in);
                % 0-injection
                dist_value_0=Vsignal_dist3;
                dist_value_p=[Vsignal_dist3(:,1),zeros(samples_window,1)];
                dist_value_n=[Vsignal_dist3(:,1),zeros(samples_window,1)];
                out=sim(program);
                 % Vector windowed extraction
                v0=out.V0pn(td1:td2,1);
                vp=out.V0pn(td1:td2,2);
                vn=out.V0pn(td1:td2,3);
                i0=out.I0pn(td1:td2,1);
                ip=out.I0pn(td1:td2,2);
                in=out.I0pn(td1:td2,3);
                % FFT of the time-windowed signals
                V0_0=fft(v0-v0_ss)/length(v0);
                Vp_0=fft(vp-vp_ss)/length(vp);
                Vn_0=fft(vn-vn_ss)/length(vn);
                I0_0=fft(i0-i0_ss)/length(i0);
                Ip_0=fft(ip-ip_ss)/length(ip);
                In_0=fft(in-in_ss)/length(in);
            for n=1:length(In_0)
                    V0pn=[V0_0(n) V0_p(n) V0_n(n)
                          Vp_0(n) Vp_p(n) Vp_n(n)
                          Vn_0(n) Vn_p(n) Vn_n(n)];
                    I0pn=[I0_0(n) I0_p(n) I0_n(n)
                          Ip_0(n) Ip_p(n) Ip_n(n)
                          In_0(n) In_p(n) In_n(n)];
                    Z_0pn=V0pn*inv(I0pn);
                    Zpp(:,n)=Z_0pn(2,2);
                    Zpn(:,n)=Z_0pn(2,3);
                    Znp(:,n)=Z_0pn(3,2);
                    Znn(:,n)=Z_0pn(3,3);
            end
            fd0=(1:samples_window)*fs;
            end
        end
        switch linear
            case 1
                if scanner_type==1
                    pn0Plot(fd0, jw1, Ym_RLC, Ya_RLC, Ypp, Ypn, Ynp, Ynn);
                    fprintf('SIaD Tool finished. Results stored in Ypp, Ypn, Ynp and Ynn.\n');
                else
                    pn0Plot(fd0, jw1, Zm_RLC, Za_RLC, Zpp, Zpn, Znp, Znn);
                    fprintf('SIaD Tool finished. Results stored in Zpp, Zpn, Znp and Znn.\n');
                end
            case 0
                if scanner_type==1
                    pn0Plot2(fd0, Ypp, Ypn, Ynp, Ynn);
                    fprintf('SIaD Tool finished. Results stored in Ypp, Ypn, Ynp and Ynn.\n');
                else
                    pn0Plot2(fd0, Zpp, Zpn, Znp, Znn);
                    fprintf('SIaD Tool finished. Results stored in Zpp, Zpn, Znp and Znn.\n');
                end
        end
end
toc
disp('---')

function out = choose(index, varargin)
    if index >= 1 && index <= numel(varargin)
        out = varargin{index};
    else
        out = 'desconocido';
    end
end


function y = multisine( frequencyLimits, fs, Ns, varargin )
%%% Multisine: a function to generate a multi-sine signal with various
%%% properties, most notably phase distribution to minimise crest-factor.

% Author: Ben Holmes, adapted from a pseudonoise signal by Maarten van
% Walstijn.
% Date: 2019/01/09
% License: All rights reserved. (until the code is cleaned up)

% Inputs
% Required:
%   - frequencyLimits: the boundaries between which all sine components
%   will fall. Each sine will fall at multiples of f0=fs/Ns, which the
%   frequency limits will be rounded to.

%   - fs: sampling frequency.

%   - Ns: signal length in samples. It is recommended for multiple periods
%   of the signal to use repmat instead of a high value of Ns as the
%   Schroeder phases are slow to calculate, and increasing Ns will increase
%   the density of sine wave components.

% Other parameters:
%   - MagnitudeResponse: the amplitude of the sine wave components. Zero
%   values should be used for all magnitudes outside of the frequency
%   limits. Default is a flat response.

%   - PhaseResponse: Either a string to select "Schroeder",
%   "NormallyDistributed", or "ZeroValued" for the different phase options,
%   or a vector of the desired phase values. Default is "Schroeder".

%   - StartAtZero: a boolean flag to indicate whether to wrap the signal
%   such that it starts at the lowest gradient zero crossing. Default true.

%   - Normalise: a boolean flag that indicates whether to normalise the
%   signal to a peak value of 1. Default true.

%   - TimeDomain: a boolean flag that indicates whether to generate the
%   signal in the time domain or frequency domain. Default false. Used for
%   debugging the IFFT.

%   - InitialPhase: a scalar element that sets the first value of the
%   Schroeder phases, ignored for other settings. Default 0.

% Output
%   y: the multi-sine output signal.

%% Input parsing

p = inputParser;

is2ElementPositiveVector =@(x) isnumeric(x) && any(size(x) == 1)...
                                            && any(size(x) == 2)...
                                            && ~any(x < 0);
                                        
addRequired(p, 'frequencyLimits', is2ElementPositiveVector);

isPositiveScalarInteger =@(x) isnumeric(x) && isscalar(x) && (round(x) == x) && x > 0;

addRequired(p, 'fs', isPositiveScalarInteger);
addRequired(p, 'Ns', isPositiveScalarInteger);

addParameter(p, 'MagnitudeResponse', false, @(x) isnumeric(x) && any(size(x) == 1) && ~any(x < 0) && sum(x)>0)

addParameter(p, 'PhaseResponse', 'Schroeder', @(x) ischar(x) || isnumeric(x))

addParameter(p, 'StartAtZero', true, @(x) islogical(x) && isscalar(x));

addParameter(p, 'Normalise', true, @(x) islogical(x) && isscalar(x));

addParameter(p, 'TimeDomain', false, @(x) islogical(x) && isscalar(x));

addParameter(p, 'InitialPhase', 0, @(x) isnumeric(x) && isscalar(x));
         
parse(p, frequencyLimits, fs, Ns, varargin{:});

%% Find frequency indices

f0 = fs/Ns;

% DC is in bin 1 so everything must start at 2
fInds = 1 + round(frequencyLimits./f0);

if any(fInds > Ns/2)
    error('Frequency limits must be beneath Nyquist.');
end

indexVector = fInds(1):fInds(2);

NN = length(indexVector);

%% Find amplitude response

if ~any(p.Results.MagnitudeResponse)
    mag = zeros(1, Ns);
    mag(indexVector) = 1./length(indexVector);
else
    if length(p.Results.MagnitudeResponse) ~= Ns
        error('Magnitude response must be the same length as the desired signal.');
    end
    
    mag = p.Results.MagnitudeResponse.^2;
    
    % Find indices at which no components should be present.
    fullIndices = (1:Ns);
    zeroValueIndices = fullIndices;
    zeroValueIndices(indexVector) = [];
    
    if any(mag(zeroValueIndices))
        warning('Non-zero magnitude values present outside of frequency limits.');
        mag(zeroValueIndices) = zeros(size(zeroValueIndices));
    end
    
    
    if sum(mag(indexVector)) ~= 1
        mag(indexVector) = mag(indexVector)./sum(mag(indexVector));
    end
end

%% Find phase response

if ischar(p.Results.PhaseResponse)
    switch p.Results.PhaseResponse
        case 'Schroeder'
            phase = schroederPhases(NN, Ns, indexVector, mag, p.Results.InitialPhase);
        case 'ZeroPhase'
            phase = zeros(1, Ns);
        case 'NormalDistribution'
            phase = randn(1, Ns);
        otherwise
            error('Phase Response string must be Schroeder, ZeroPhase, or NormalDistribution.');
    end
else
    phase = p.Results.PhaseResponse;
    if ~any(size(phase) == 1) || ~any(size(phase) == Ns)
        error('Phase response must be Ns x 1.');
    end
end


%% Generate multisine signal

% Switch between time and frequency domain generation
if p.Results.TimeDomain
    y = zeros(1, Ns);
    t = (0:Ns-1)./fs;
    for nn=1:NN
        mm = indexVector(nn);
        y = y + sqrt(mag(mm)/2)*cos(2*pi*f0*(mm-1)*t + phase(mm));
    end
else
    % Frequency domain representation
    Y = sqrt(mag/2).*exp(sqrt(-1).*phase);

    % IFFT to time domain
    y = ifft(forceFFTSymmetry(Y))*(Ns/2);
end

%% Normalise peak abs value to 1

if p.Results.Normalise
    y = y./max(abs(y));
end

%% Find zero crossing closest to zero 

% Heuristic method of finding minimum gradient zero crossing
if p.Results.StartAtZero
    ySign = y>0;

    zeroInds = find((ySign(2:end) ~= ySign(1:end-1)));

    % Find the index with the smallest gradient around the zero crossing.
    zeroGrad = zeros(1,length(zeroInds));
    for nn=1:length(zeroInds)
        zeroGrad(nn) = abs(y(zeroInds(nn)) - y(zeroInds(nn)+1));
    end
    [~, minInd] = min(zeroGrad);

    yWrapped = [y(zeroInds(minInd):end) y(1:zeroInds(minInd)-1)];

    y = yWrapped;
end

end

function phase = schroederPhases(nComponents, Ns, indexVector, magnitude, p1)
% Generate phases as defined in "Synthesis of Low-Peak-Factor Signals and
% Binary Sequences With Low Autocorrelation" by M. R. Schroeder

% Preallocate vector
phase = zeros(1, Ns);

% Bin 1 is DC, so bin 2 is the phase of the first component
phase(2) = p1;

% Iterate over phase values using Schroeder's algorithm
for nn=2:nComponents
    ll=1:(nn-1);
    phase(indexVector(nn)) = phase(2) -2*pi*sum((nn-ll).*magnitude(indexVector(ll)));
end

end

function Y = forceFFTSymmetry(X)
% forceFFTSymmetry  A function to force conjugate symmetry on an FFT such that when an
% IFFT is performed the result is a real signal.

% The function has been written to replace MATLAB's ifft(X,'symmetric'), as this function
% is not compatible with MATLAB Coder.
  
Y = X;
XStartFlipped = fliplr(X(2:floor(end/2)));
Y(ceil(end/2)+2:end) = real(XStartFlipped) - sqrt(complex(-1))*imag(XStartFlipped);

end

function ABCPlot(fd0, Yabc, Yphases, jw1)
    % ABCPlot: Generates a 6x3 frequency response plot.
    %
    % Input parameters:
    % - fd0: Vector of scanned frequencies (Hz)
    % - Yabc: 3x3xn matrix with scanned responses
    % - Yphases: 3x3xm matrix with theoretical responses
    % - jw1: Vector of theoretical frequencies in the complex domain

    % Define x-axis limits
    low_axis = fd0(1);
    up_axis = fd0(end);

    % Global plot settings
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Create figure
    figure;

    % Iterate over row and column combinations (3x3)
    for i = 1:3
        for j = 1:3
            % Extract scanned and theoretical responses for current position
            Ya_scan = squeeze(Yabc(i, j, :));       % Scanned response (vector)
            Yphase_theo = squeeze(Yphases(i, j, :)); % Theoretical response (vector)

            % Subplot indices
            mag_idx = 2*(i-1) * 3 + j;   % Magnitude index
            phase_idx = mag_idx + 3;    % Phase index

            % Subplot for magnitude
            subplot(6, 3, mag_idx);
            semilogx(imag(jw1)/(2*pi), 20*log10(abs(Yphase_theo)), 'k'); % Theoretical response
            hold on;
            semilogx(fd0, 20*log10(abs(Ya_scan)), 'rx'); % Measured response
%             title(['Y_{', char(64+i), char(64+j), '}(s) Magnitude']);
            ylabel('Magnitude (dB)');
            xlim([low_axis up_axis]);
            grid on; grid minor;

            % Subplot for phase
            subplot(6, 3, phase_idx);
            semilogx(imag(jw1)/(2*pi), (180/pi)*angle(Yphase_theo), 'k'); % Theoretical response
            hold on;
            semilogx(fd0, (180/pi)*angle(Ya_scan), 'rx'); % Measured response
%             title(['Y_{', char(64+i), char(64+j), '}(s) Phase']);
            ylabel('Phase (deg)');
            xlabel('Frequency (Hz)');
            xlim([low_axis up_axis]);
            grid on; grid minor;
        end
    end
end

function ABCPlot2(fd0, Yabc)
    % ABCPlot2: Generates a 6x3 plot of scanned frequency responses.
    %
    % Input parameters:
    % - fd0: Vector of scanned frequencies (Hz)
    % - Yabc: 3x3xn matrix with scanned responses

    % Define x-axis limits
    low_axis = fd0(1);
    up_axis = fd0(end);

    % Global plot settings
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Create figure
    figure;

    % Iterate over row and column combinations (3x3)
    for i = 1:3
        for j = 1:3
            % Extract scanned response for current position
            Ya_scan = squeeze(Yabc(i, j, :)); % Scanned response (vector)
            
            % Ensure Ya_scan is a column vector
            Ya_scan = Ya_scan(:);

            % Subplot indices
            mag_idx = 2*(i-1) * 3 + j;   % Magnitude index
            phase_idx = mag_idx + 3;    % Phase index

            % Subplot for magnitude
            subplot(6, 3, mag_idx);
            semilogx(fd0, 20*log10(abs(Ya_scan)), 'r-'); % Measured response
            title(['Y_{', char(64+i), char(64+j), '}(s) Magnitude']);
            ylabel('Magnitude (dB)');
            xlim([low_axis up_axis]);
            grid on; grid minor;

            % Subplot for phase
            subplot(6, 3, phase_idx);
            semilogx(fd0, (180/pi)*angle(Ya_scan), 'r-'); % Measured response
            title(['Y_{', char(64+i), char(64+j), '}(s) Phase']);
            ylabel('Phase (deg)');
            xlabel('Frequency (Hz)');
            xlim([low_axis up_axis]);
            grid on; grid minor;
        end
    end

    % Overall title for all subplots
    sgtitle('Frequency Response Plots (Scanned Data)');
end

function qd0Plot(fd0, jw1, Ym_Th, Ya_Th, Yqq, Yqd, Ydq, Ydd)
    % Limits of x axis
    low_axis = fd0(1);
    up_axis = fd0(end);

    % Global configurations
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Create figure
    figure;

    % Subplot 1: Magnitude of Yqq
    subplot(4, 2, 1);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(1, 1, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Yqq)), 'rx');
    title('Yqq(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 2: Phase of Yqq
    subplot(4, 2, 3);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(1, 1, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Yqq), 'rx');
    ylabel('Phase (deg)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 3: Magnitude de Yqd
    subplot(4, 2, 2);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(1, 2, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Yqd)), 'rx');
    legend({'Linear: state space model', 'dq0: Frequency scan'}, 'Location', 'southwest', 'Orientation', 'vertical');
    title('Yqd(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 4: Phase de Yqd
    subplot(4, 2, 4);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(1, 2, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Yqd), 'rx');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 5: Magnitude de Ydq
    subplot(4, 2, 5);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(2, 1, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Ydq)), 'rx');
    title('Ydq(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 6: Phase de Ydq
    subplot(4, 2, 7);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(2, 1, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Ydq), 'rx');
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 7: Magnitude de Ydd
    subplot(4, 2, 6);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(2, 2, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Ydd)), 'rx');
    title('Ydd(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 8: Phase de Ydd
    subplot(4, 2, 8);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(2, 2, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Ydd), 'rx');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Title
    sgtitle('Frequency Response Plots');
end

function qd0Plot2(fd0, Yqq, Yqd, Ydq, Ydd)
    % Definir límites del eje x
    low_axis = fd0(1);
    up_axis = fd0(end);

    % Configuraciones globales para las gráficas
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Crear la figura
    figure;

    % Subplot 1: Magnitud de Yqq
    subplot(4, 2, 1);
    semilogx(fd0, 20*log10(abs(Yqq)), 'r-'); % Respuesta medida
    title('Yqq(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 2: Fase de Yqq
    subplot(4, 2, 3);
    semilogx(fd0, (180/pi) * angle(Yqq), 'r-'); % Respuesta medida
    ylabel('Phase (deg)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 3: Magnitud de Yqd
    subplot(4, 2, 2);
    semilogx(fd0, 20*log10(abs(Yqd)), 'r-'); % Respuesta medida
    legend({'qd frequency scan'}, 'Location', 'southwest', 'Orientation', 'vertical');
    title('Yqd(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 4: Fase de Yqd
    subplot(4, 2, 4);
    semilogx(fd0, (180/pi) * angle(Yqd), 'r-'); % Respuesta medida
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 5: Magnitud de Ydq
    subplot(4, 2, 5);
    semilogx(fd0, 20*log10(abs(Ydq)), 'r-'); % Respuesta medida
    title('Ydq(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 6: Fase de Ydq
    subplot(4, 2, 7);
    semilogx(fd0, (180/pi) * angle(Ydq), 'r-'); % Respuesta medida
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 7: Magnitud de Ydd
    subplot(4, 2, 6);
    semilogx(fd0, 20*log10(abs(Ydd)), 'r-'); % Respuesta medida
    title('Ydd(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 8: Fase de Ydd
    subplot(4, 2, 8);
    semilogx(fd0, (180/pi) * angle(Ydd), 'r-'); % Respuesta medida
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Título general
    sgtitle('Frequency Response Plots');
end

function pn0Plot(fd0, jw1, Ym_Th, Ya_Th, Ypp, Ypn, Ynp, Ynn)
    % Definir límites del eje x
    low_axis = fd0(1);
    up_axis = fd0(end);

    % Configuraciones globales para las gráficas
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Crear la figura
    figure;

    % Subplot 1: Magnitud de Ypp
    subplot(4, 2, 1);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(1, 1, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Ypp)), 'rx');
    title('Ypp(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 2: Fase de Ypp
    subplot(4, 2, 3);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(1, 1, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Ypp), 'rx');
    ylabel('Phase (deg)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 3: Magnitud de Ypn
    subplot(4, 2, 2);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(1, 2, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Ypn)), 'rx');
    legend({'Linear: state space model', 'pn0: Frequency scan'}, 'Location', 'southwest', 'Orientation', 'vertical');
    title('Ypn(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 4: Fase de Ypn
    subplot(4, 2, 4);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(1, 2, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Ypn), 'rx');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 5: Magnitud de Ynp
    subplot(4, 2, 5);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(2, 1, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Ynp)), 'rx');
    title('Ynp(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 6: Fase de Ynp
    subplot(4, 2, 7);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(2, 1, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Ynp), 'rx');
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 7: Magnitud de Ynn
    subplot(4, 2, 6);
    semilogx(imag(jw1)/(2*pi), 20*log10(squeeze(Ym_Th(2, 2, :))), 'k');
    hold on;
    semilogx(fd0, 20*log10(abs(Ynn)), 'rx');
    title('Ynn(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 8: Fase de Ynn
    subplot(4, 2, 8);
    semilogx(imag(jw1)/(2*pi), squeeze(Ya_Th(2, 2, :)), 'k');
    hold on;
    semilogx(fd0, (180/pi) * angle(Ynn), 'rx');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Título general
    sgtitle('Frequency Response Plots');
end

function pn0Plot2(fd0, Ypp, Ypn, Ynp, Ynn)
    % Definir límites del eje x
    low_axis = fd0(1);
    up_axis = fd0(end);

    % Configuraciones globales para las gráficas
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Crear la figura
    figure;

    % Subplot 1: Magnitud de Ypp
    subplot(4, 2, 1);
    semilogx(fd0, 20*log10(abs(Ypp)), 'r-'); % Respuesta medida
    title('Ypp(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 2: Fase de Ypp
    subplot(4, 2, 3);
    semilogx(fd0, (180/pi) * angle(Ypp), 'r-'); % Respuesta medida
    ylabel('Phase (deg)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 3: Magnitud de Ypn
    subplot(4, 2, 2);
    semilogx(fd0, 20*log10(abs(Ypn)), 'r-'); % Respuesta medida
    legend({'pn0: Frequency scan'}, 'Location', 'southwest', 'Orientation', 'vertical');
    title('Ypn(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 4: Fase de Ypn
    subplot(4, 2, 4);
    semilogx(fd0, (180/pi) * angle(Ypn), 'r-'); % Respuesta medida
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 5: Magnitud de Ynp
    subplot(4, 2, 5);
    semilogx(fd0, 20*log10(abs(Ynp)), 'r-'); % Respuesta medida
    title('Ynp(s)');
    ylabel('Magnitude (dB)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 6: Fase de Ynp
    subplot(4, 2, 7);
    semilogx(fd0, (180/pi) * angle(Ynp), 'r-'); % Respuesta medida
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 7: Magnitud de Ynn
    subplot(4, 2, 6);
    semilogx(fd0, 20*log10(abs(Ynn)), 'r-'); % Respuesta medida
    title('Ynn(s)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Subplot 8: Fase de Ynn
    subplot(4, 2, 8);
    semilogx(fd0, (180/pi) * angle(Ynn), 'r-'); % Respuesta medida
    xlabel('Frequency (Hz)');
    xlim([low_axis up_axis]);
    grid on; grid minor;

    % Título general
    sgtitle('Frequency Response Plots');
end
