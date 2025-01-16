function [Zm_RLC, Za_RLC] = Zmatrix_qd_ss(samples, Ym_RLC, Ya_RLC)
    % FUNCION QUE COMPUTA MATRICES Zm_RLC Y Za_RLC
    % Argumentos de entrada:
    % samples    - Número de muestras
    % Ym_RLC     - Tensor 3D con matrices Ym
    % Ya_Th_gfl  - Tensor 3D con matrices Ya (ángulos en grados)
    %
    % Salidas:
    % Zm_RLC     - Tensor 3D con magnitudes de impedancia
    % Za_RLC     - Tensor 3D con ángulos de impedancia (en grados)
    
    % Inicialización de tensores de salida
    Zm_RLC = zeros(2, 2, samples);
    Za_RLC = zeros(2, 2, samples);

    % Loop principal
    for n = 1:samples
        % Extraer matrices Ym y Ya
        Ym = squeeze(Ym_RLC(1:2, 1:2, n));
        Ya = squeeze(Ya_RLC(1:2, 1:2, n));
        
        % Calcular elementos de la matriz Y modificada
        Y11 = Ym(1, 1) * (cosd(Ya(1, 1)) + 1i * sind(Ya(1, 1)));
        Y12 = Ym(1, 2) * (cosd(Ya(1, 2)) + 1i * sind(Ya(1, 2)));
        Y21 = Ym(2, 1) * (cosd(Ya(2, 1)) + 1i * sind(Ya(2, 1)));
        Y22 = Ym(2, 2) * (cosd(Ya(2, 2)) + 1i * sind(Ya(2, 2)));
        
        % Construir matriz Y0 y calcular su inversa Z0
        Z0 = inv([Y11, Y12; Y21, Y22]);
        
        % Obtener magnitudes
        Z11m = abs(Z0(1, 1));
        Z12m = abs(Z0(1, 2));
        Z21m = abs(Z0(2, 1));
        Z22m = abs(Z0(2, 2));
        
        % Obtener ángulos en grados
        Z11a = rad2deg(angle(Z0(1, 1)));
        Z12a = rad2deg(angle(Z0(1, 2)));
        Z21a = rad2deg(angle(Z0(2, 1)));
        Z22a = rad2deg(angle(Z0(2, 2)));
        
        % Almacenar resultados en tensores
        Zm_RLC(:, :, n) = [Z11m, Z12m; Z21m, Z22m];
        Za_RLC(:, :, n) = [Z11a, Z12a; Z21a, Z22a];
    end
end
