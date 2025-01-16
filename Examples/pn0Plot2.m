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
