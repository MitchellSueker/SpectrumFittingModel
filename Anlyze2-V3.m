clear
global cnt
cnt = 0;

disp('analyze skin1')

%lookskin1  % loads: nm, skin, std, R
newlookskin
nm0 = nm;  % store original wavelengths

numSpectra = size(data, 2);  % number of spectra to fit

% Load absorption coefficients
loadTakataniGraham     % nmTG, muaoxyTG, muadeoxyTG
% loadwaterHaleQuerry    % waterHaleQuerry
loadwaterS
nmHQ = waterHaleQuerry(:,1);
muaHQ = waterHaleQuerry(:,2);
clear waterHaleQuerry
fprintf('\n[DEBUG] Measured wavelength range: %.1f – %.1f nm\n', min(nm), max(nm));
fprintf('[DEBUG] Hemoglobin range (nmTG): %.1f – %.1f nm\n', min(nmTG), max(nmTG));
fprintf('[DEBUG] Water range (nmHQ): %.1f – %.1f nm\n\n', min(nmHQ), max(nmHQ));
drawnow

% Get wavelength indices for fitting
% ✅ Define safe fitting range (450–1700 nm) to avoid interpolation NaNs
% --- Find index where nm >= 450 ---
i = 1;
while i <= length(nm) && nm(i) < 450
    i = i + 1;
end
i450 = i;

% --- Find index where nm >= 1700 ---
i = 1;
while i <= length(nm) && nm(i) < 1700
    i = i + 1;
end
i1700 = i - 1;  % stop before crossing 1700 or array end

% --- Define fitting range ---
u = (i450:i1700);  % valid wavelengths only (450–1700 nm)


% Emphasize 500–600 nm region for hemoglobin fitting
i = 1;
while nm(u(i)) < 500; i = i + 1; end; u500 = i;
while nm(u(i)) < 600; i = i + 1; end; u600 = i;
v = (u500:u600);  % Hb-sensitive region (subset of u)


% Normalize reflectance
M = R;  % M = Mskin / Mstd

% ----- Optimization Setup -----

start = [
    0.3, 0.01, 0.9, 5, 0.02, 1.0, 0.1, 0.006;
    0.4, 0.015, 0.95, 8, 0.03, 1.0, 0.1, 0.006;
    0.5, 0.02, 0.85, 10, 0.04, 1.0, 0.1, 0.006;
    0.6, 0.005, 0.92, 12, 0.05, 1.0, 0.1, 0.006;
    0.35, 0.008, 0.88, 6, 0.03, 1.0, 0.1, 0.006;
];  % [W B S a Mel const nmOff Lepi]
lb =    [0.2,   0.0003, 0.2,  0.1, 0.001, 0.5, 0, 0.003];
ub =    [0.7,   0.05, 1 , 30, 0.5,   2.0, 0,  0.01];

options = optimoptions('lsqnonlin', ...
    'Display','iter', ...
    'MaxIterations',20000, ...
    'TolX',1e-10, ...
    'TolFun',1e-10);

% ----- Run Optimization -----
[result, resnorm, residual] = lsqnonlin(@(x) fitspectrum2_residual(x, nm(u), M(u), ...
    nmTG, muaoxyTG, muadeoxyTG, nmHQ, muaHQ, u, v), ...
    start, lb, ub, options);

% ----- Check if optimizer updated anything -----
disp('Optimized parameters:')
disp(result)

if all(abs(result - start) < 1e-6)
    disp('⚠️ WARNING: Optimization did not update parameters!')
end

% ----- Extract and use final parameters -----
W     = result(1);
B     = result(2);
S     = result(3);
a     = result(4);
Mel   = result(5);
const = result(6);
nmOff = result(7);
Lepi  = result(8);

nm = nm0 + nmOff;

% ----- Recompute model using final result -----
muaoxy = interp1(nmTG, muaoxyTG, nm);
muadeoxy = interp1(nmTG, muadeoxyTG, nm);
muawater = interp1(nmHQ, muaHQ, nm);
% muaoxy = interp1(nmTG, muaoxyTG, nm, 'pchip');
% muadeoxy = interp1(nmTG, muadeoxyTG, nm, 'pchip');
% muawater = interp1(nmHQ, muaHQ, nm, 'pchip');
muamel = 6.6e11 * nm.^-3.33;

Mie = 4.59e3 * nm.^-0.913;
Ray = 1.74e12 * nm.^-4;
musp = a * (Mie + Ray);

n = 1.4;
Tepi = exp(-Mel * muamel * 2 * Lepi);
mua  = B * (S*muaoxy + (1-S)*muadeoxy) + W * muawater;
pRwMB = const * Tepi .* getRdFarrell(mua, musp, n);

% ----- Plotting -----
figure(3); clf
set(figure(3), 'position', [645, 40, 577, 757], 'color', 'w')
sz = 18;

% Top: Measured vs Modeled
subplot(2,1,1); hold off
plot(nm(u), M(u), 'ko','linewidth',1)
hold on
plot(nm, pRwMB, 'r-','linewidth',2)
set(gca,'fontsize',sz,'linewidth',2)
xlabel('wavelength [nm]')
ylabel('Reflectance')

title('Measured vs. Modeled Spectrum')
axis([450 1700 0 1])
x = 510; ymax = 1; dy = .07;
text(x, ymax - dy, sprintf('W = %0.1f', W*100),'fontsize',sz)
text(x, ymax - 2*dy, sprintf('B = %0.4f', B),'fontsize',sz)
text(x, ymax - 3*dy, sprintf('S = %0.3f', S),'fontsize',sz)
text(x, ymax - 4*dy, sprintf('a = %0.3f', a),'fontsize',sz)
text(x, ymax - 5*dy, sprintf('const = %0.2f', const),'fontsize',sz)
text(x, ymax - 6*dy, sprintf('Mel = %0.3f', Mel),'fontsize',sz)
text(x, ymax - 7*dy, sprintf('Lepi = %0.3f', Lepi),'fontsize',sz)

% Bottom: absorption components
% subplot(2,1,2); hold off
% semilogy(nm, W*muawater, 'b--','linewidth',2)
% hold on
% plot(nm, B*S*muaoxy, 'r-','linewidth',2)
% plot(nm, B*(1-S)*muadeoxy, 'b-','linewidth',2)
% plot(nm, Mel*muamel, 'k-','linewidth',2)
% plot(nm, mua,'m-','linewidth',2)
% plot(nm, Tepi,'k-','linewidth',2)
% set(gca,'fontsize',sz,'linewidth',2)
% xlabel('wavelength [nm]')
% ylabel('\mu_a [cm^-^1]')
% axis([450 1700 1e-3 1])