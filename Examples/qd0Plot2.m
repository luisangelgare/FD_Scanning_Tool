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
