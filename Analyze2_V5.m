clear
global cnt
cnt = 0;

disp('analyze skin1')

newlookskin
nm0 = nm;  % store original wavelengths
numSpectra = size(data, 2);  % number of spectra to fit
W_values = [];

for i = 2:numSpectra
    R = data{:, i};  % get the i-th reflectance spectrum
    window_size = 21;  % Must be odd
    poly_order = 1;  % First-order polynomial smoothing
    R = sgolayfilt(R, poly_order, window_size);

    % Load absorption coefficients
    loadTakataniGraham
    loadwaterS
    nmHQ = waterHaleQuerry(:,1);
    muaHQ = waterHaleQuerry(:,2);

    % Define fitting range
    idx = 1;
    while idx <= length(nm) && nm(idx) < 450
        idx = idx + 1;
    end
    i450 = idx;

    idx = 1;
    while idx <= length(nm) && nm(idx) < 1700
        idx = idx + 1;
    end
    i1700 = idx - 1;
    u = (i450:i1700);

    idx = 1;
    while nm(u(idx)) < 500; idx = idx + 1; end; u500 = idx;
    while nm(u(idx)) < 600; idx = idx + 1; end; u600 = idx;
    v = (u500:u600);

    M = R;

    % Initial guess set
    start = [
        0.3, 0.01, 0.9, 5, 0.02, 1.0, 0.1, 0.006;
        0.4, 0.015, 0.95, 8, 0.03, 1.0, 0.1, 0.006;
        0.5, 0.02, 0.85, 10, 0.04, 1.0, 0.1, 0.006;
        0.6, 0.005, 0.92, 12, 0.05, 1.0, 0.1, 0.006;
        0.35, 0.008, 0.88, 6, 0.03, 1.0, 0.1, 0.006;
    ];

    lb = [0.2, 0.0003, 0.2, 0.1, 0.001, 0.5, 0, 0.003];
    ub = [0.7, 0.05, 1, 30, 0.5, 2.0, 0, 0.01];

    options = optimoptions('lsqnonlin', ...
        'Display', 'off', ...
        'MaxIterations', 20000, ...
        'TolX', 1e-2, ...
        'TolFun', 1e-2);

    % Define weight ranges
    weight_Hb_all = 0:0.5:3;
    weight_H2O_all = 0:0.5:3;

    bestW = NaN;
    minResnorm = Inf;
    bestResult = NaN(1, 8);
    bestWeights = NaN(1, 2);

    % Loop through weight combinations
    for wHb = weight_Hb_all
        for wH2O = weight_H2O_all
            weights_vector = [wHb, wH2O];

            % Try multiple starting points
            for s = 1:size(start, 1)
                [result, resnorm, ~] = lsqnonlin(@(x) fitspectrum2_residual_loopweightspectra(x, nm(u), M(u), ...
                    nmTG, muaoxyTG, muadeoxyTG, nmHQ, muaHQ, u, v, weights_vector), ...
                    start(s, :), lb, ub, options);

                if resnorm < minResnorm
                    minResnorm = resnorm;
                    bestResult = result;
                    bestW = result(1);
                    bestWeights = weights_vector;
                end
            end
        end
    end

    % Save best W
    W_values = [W_values, bestW];

    % Optional: Display result for this spectrum
    fprintf('Spectrum %d: Best W = %.4f | Hb weight = %.2f | H2O weight = %.2f | Resnorm = %.4f\n', ...
        i - 1, bestW, bestWeights(1), bestWeights(2), minResnorm);
end
