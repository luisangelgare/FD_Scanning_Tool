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
scanner_selector=2;

% 1 -> linear response exists
% 0 -> no linear response is available
linear=0;

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
 
%% Comprehensive frequency scanning tool methodology

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