function ABCPlot(fd0, Yabc, Yphases, jw1)
    % ABCPlot: Genera una gráfica de 6x3 de respuestas de frecuencia.
    %
    % Parámetros de entrada:
    % - fd0: Vector de frecuencias escaneadas (Hz)
    % - Yabc: Matriz 3x3xn con respuestas escaneadas
    % - Yphases: Matriz 3x3xm con respuestas teóricas
    % - jw1: Vector de frecuencias teóricas en el dominio complejo

    % Definir límites del eje x
    low_axis = fd0(1);
    up_axis = fd0(end);

    % Configuraciones globales para las gráficas
    set(0, 'defaultAxesFontSize', 14);
    set(0, 'DefaultLineLineWidth', 1.5);

    % Crear la figura
    figure;

    % Iterar sobre las combinaciones de filas y columnas (3x3)
    for i = 1:3
        for j = 1:3
            % Obtener la respuesta escaneada y teórica para la posición actual
            Ya_scan = squeeze(Yabc(i, j, :)); % Respuesta escaneada (vector)
            Yphase_theo = squeeze(Yphases(i, j, :)); % Respuesta teórica (vector)

            % Índices de los subplots
            mag_idx = 2*(i-1) * 3 + j; % Índice de la magnitud
            phase_idx = mag_idx + 3; % Índice de la fase

            % Subplot para la magnitud
            subplot(6, 3, mag_idx);
            semilogx(imag(jw1)/(2*pi), 20*log10(abs(Yphase_theo)), 'k'); % Respuesta teórica
            hold on;
            semilogx(fd0, 20*log10(abs(Ya_scan)), 'rx'); % Respuesta medida
%             title(['Y_{', char(64+i), char(64+j), '}(s) Magnitude']);
            ylabel('Magnitude (dB)');
            xlim([low_axis up_axis]);
            grid on; grid minor;

            % Subplot para la fase
            subplot(6, 3, phase_idx);
            semilogx(imag(jw1)/(2*pi), (180/pi)*angle(Yphase_theo), 'k'); % Respuesta teórica
            hold on;
            semilogx(fd0, (180/pi)*angle(Ya_scan), 'rx'); % Respuesta medida
%             title(['Y_{', char(64+i), char(64+j), '}(s) Phase']);
            ylabel('Phase (deg)');
            xlabel('Frequency (Hz)');
            xlim([low_axis up_axis]);
            grid on; grid minor;
        end
    end

    % Título general para todas las gráficas
%     sgtitle('Frequency Response Plots');
end


% function ABCPlot(fd0, Ya, Yb, Yc, Yphase, jw1)
%     % plot_frequency_response: Genera gráficas de respuesta de frecuencia.
%     %
%     % Parámetros de entrada:
%     % - fd0: Vector de frecuencias para escanear (Hz)
%     % - Ya: Magnitud de Yaa
%     % - Yb: Magnitud de Ybb
%     % - Yc: Magnitud de Ycc
%     % - Yphase: Fase para la respuesta teórica
%     % - jw1: Vector para las frecuencias teóricas en el dominio complejo
%     
%     % Definir límites del eje x
%     low_axis = fd0(1);
%     up_axis = fd0(end);
% 
%     % Configuraciones globales para las gráficas
%     set(0, 'defaultAxesFontSize', 14);
%     set(0, 'DefaultLineLineWidth', 1.5);
% 
%     % Crear la figura
%     figure;
% 
%     % Subplot 1: Magnitud de Yaa
%     subplot(2, 3, 1);
%     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Yphase))), 'k'); % Respuesta teórica
%     hold on;
%     semilogx(fd0, 20*log10(abs(Ya)), 'rx'); % Respuesta medida
%     title('Y_{AA}(s)');
%     ylabel('Magnitude (dB)');
%     xlim([low_axis up_axis]);
%     grid on; grid minor;
% 
%     % Subplot 2: Fase de Yaa
%     subplot(2, 3, 4);
%     semilogx(imag(jw1)/(2*pi), (180/pi)*angle(squeeze(Yphase)), 'k'); % Respuesta teórica
%     hold on;
%     semilogx(fd0, (180/pi) * angle(Ya), 'rx'); % Respuesta medida
%     ylabel('Phase (deg)');
%     xlabel('Frequency (Hz)');
%     xlim([low_axis up_axis]);
%     grid on; grid minor;
% 
%     % Subplot 3: Magnitud de Ybb
%     subplot(2, 3, 2);
%     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Yphase))), 'k'); % Respuesta teórica
%     hold on;
%     semilogx(fd0, 20*log10(abs(Yb)), 'rx'); % Respuesta medida
%     title('Y_{BB}(s)');
%     xlim([low_axis up_axis]);
%     grid on; grid minor;
% 
%     % Subplot 4: Fase de Ybb
%     subplot(2, 3, 5);
%     semilogx(imag(jw1)/(2*pi), (180/pi)*angle(squeeze(Yphase)), 'k'); % Respuesta teórica
%     hold on;
%     semilogx(fd0, (180/pi) * angle(Yb), 'rx'); % Respuesta medida
%     xlabel('Frequency (Hz)');
%     xlim([low_axis up_axis]);
%     grid on; grid minor;
% 
%     % Subplot 5: Magnitud de Ycc
%     subplot(2, 3, 3);
%     semilogx(imag(jw1)/(2*pi), 20*log10(abs(squeeze(Yphase))), 'k'); % Respuesta teórica
%     hold on;
%     semilogx(fd0, 20*log10(abs(Yc)), 'rx'); % Respuesta medida
%     title('Y_{CC}(s)');
%     legend({'Theoretical response', 'ABC frequency scan'}, 'Location', 'southwest', 'Orientation', 'vertical');
%     xlim([low_axis up_axis]);
%     grid on; grid minor;
% 
%     % Subplot 6: Fase de Ycc
%     subplot(2, 3, 6);
%     semilogx(imag(jw1)/(2*pi), (180/pi)*angle(squeeze(Yphase)), 'k'); % Respuesta teórica
%     hold on;
%     semilogx(fd0, (180/pi) * angle(Yc), 'rx'); % Respuesta medida
%     xlabel('Frequency (Hz)');
%     xlim([low_axis up_axis]);
%     grid on; grid minor;
% 
%     % Título general
%     sgtitle('Frequency Response Plots');
% end
