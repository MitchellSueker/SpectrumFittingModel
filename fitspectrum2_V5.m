function res = fitspectrum2_residual(x, nm, M, nmTG, muaoxyTG, muadeoxyTG, nmHQ, muaHQ, u, v, weights_vector)

% Unpack parameters
W     = x(1);
B     = x(2);
S     = x(3);
a     = x(4);
Mel   = x(5);
const = x(6);
nmOff = x(7);
Lepi  = x(8);

% Shift wavelengths
nm_shifted = nm + nmOff;

% Interpolate absorption coefficients
muaoxy    = interp1(nmTG, muaoxyTG, nm, 'pchip');
muadeoxy  = interp1(nmTG, muadeoxyTG, nm, 'pchip');
muawater  = interp1(nmHQ, muaHQ, nm, 'pchip');
muamel    = 6.6e11 * nm_shifted.^-3.33;

% Safety check
if any(isnan(muaoxy)) || any(isnan(muadeoxy)) || any(isnan(muawater))
    disp('❌ Interpolation failed: NaN in chromophore data');
    res = ones(size(M)) * 1e6;
    return
end

% Scattering (Mie + Rayleigh)
Mie = 4.59e3 * nm_shifted.^-0.913;
Ray = 1.74e12 * nm_shifted.^-4;
musp = a * (Mie + Ray);

% Absorption and transmission
mua  = B * (S * muaoxy + (1 - S) * muadeoxy) + W * muawater;
Tepi = exp(-Mel * muamel * 2 * Lepi);

if any(isnan(mua)) || any(isnan(musp)) || any(isnan(Tepi))
    disp('❌ NaN in mua, musp, or Tepi');
    res = ones(size(M)) * 1e6;
    return
end

% Compute model reflectance
try
    Rmodel = const * Tepi .* getRdFarrell(mua, musp, 1.4);
catch
    disp('❌ getRdFarrell failed');
    res = ones(size(M)) * 1e6;
    return
end

% Check model validity
if all(isnan(Rmodel)) || all(Rmodel == Rmodel(1))
    disp('❌ Model output is constant or NaN');
    res = ones(size(M)) * 1e10;
    return
end

% Compute residuals
residual_full = M - Rmodel;

% Apply custom weights
weights = ones(size(nm));  % default = 1
weights(nm >= 450 & nm <= 600) = weights_vector(1);     % Hb weight
weights(nm >= 1400 & nm <= 1500) = weights_vector(2);   % Water weight

% Apply weights only in fitting region
res = weights(u) .* residual_full(u);

% Penalize unrealistic S values
if x(3) > 0.98
    penalty = 10 * (x(3) - 0.98);
    res = res + penalty;
end

% Optional debug plot every 50 calls
persistent callCount
if isempty(callCount)
    callCount = 1;
else
    callCount = callCount + 1;
end

if callCount <= 10 || mod(callCount, 50) == 0
    fprintf('[DEBUG] Call #%d\n', callCount);
    fprintf('x = %s\n', num2str(x, '%.4f '));
    fprintf('Resnorm = %.4f\n\n', norm(res)^2);

    figure(99); clf;
    plot(nm(u), M(u), 'ko'); hold on;
    plot(nm(u), Rmodel(u), 'r-', 'linewidth', 2);
    title(sprintf('Iteration %d: Measured vs. Modeled', callCount));
    xlabel('Wavelength [nm]');
    ylabel('Reflectance');
    legend('Measured', 'Model');
    axis([450 1700 0 1]);
    drawnow;
end

end
