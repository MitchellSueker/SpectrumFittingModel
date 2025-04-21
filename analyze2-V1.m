% analyze2.m
% Uses
%	lookskin1.m
%       Reference and White table data
%	loadTakataniGraham.m
%	loadwaterHaleQuerry.m
%	fitspectrum2.m
%	getRdfarrel.m

clear
global cnt
cnt = 0;

disp('analyze skin1')


lookskin1  % loads spectrum --> nm, skin, std, R
%Mitch: Automatically get data from excel sheet

nm0 = nm;

%%%%%%
% Load LIBRARY of spectra
%%%%
loadTakataniGraham  % hemoglobin (oxy, deoxy)
% Mitch: check recent oxy; deoxy Hb, Melanin and water absorption
% coefficients
%
loadwaterHaleQuerry  % water
nmHQ = waterHaleQuerry(:,1);
muaHQ = waterHaleQuerry(:,2)
clear waterHaleQuerry
%%%%%%%
drawnow

r = 0.20; % cm = separation of detector and illumination fibers not used in PIPA

i = 1; while nm(i)<450; i = i+1; end; i450 = i;
i = 1; while nm(i)<1700; i = i+1; end; i1700 = i;
u = (i450:i1700);
i=1; while nm(u(i))<500, i=i+1;end; u500 = i;while nm(u(i))<600, i=i+1;end; u600 = i;
v = (u500:u600); % for use as nm(u(v)) in emphasizing fit to the 500-600 nm range, to specify oxy/deoxy hemoglobin


%%%%%%
% fit data
%%%%
M = R;  % repeating what was done in lookskin1.m

W = 0.30; % water content

B      = 0.002; % mean blood volume fraction in skin
S      = 0.95; % tissue oxygen saturation (mixture of arterial and venous blood)
a      = 1.0;   % scale the skin scattering (Rayleigh + Mie)
Mel    = 0.1;  % aveage volume fraction of melanosomes in 60-um-thick epidermis
const  = 0.98;  % const = Rstd*Gstd/Gskin
nmOff  = 0;     % correct for any offset error in spectrometer wavelength registration
Lepi=0.0060; % epidermis thickness
start  = [W B S a Mel const nmOff, Lepi];

% It sometimes helps to run the fitting several times, and perturb the
% blood content a little to force a readjustment.
% Mitch: check a better optimization models for a broader search/different
% intial gusses
options = optimset('MaxIter', 100000, 'TolX', 0.1); 
start = fminsearch(@(x) fitspectrum2(x, nm(u), M(u), nmTG, muaoxyTG, muadeoxyTG, nmHQ, muaHQ, u, v), start, options);
start = fminsearch(@(x) fitspectrum2(x, nm(u), M(u), nmTG, muaoxyTG, muadeoxyTG, nmHQ, muaHQ, u, v), start, options);
start(1) = start(1)/2;
result = fminsearch(@(x) fitspectrum2(x, nm(u), M(u), nmTG, muaoxyTG, muadeoxyTG, nmHQ, muaHQ, u, v), start, options);
W     = result(1);
B     = result(2);
S     = result(3);
a     = result(4);
Mel   = result(5);
const = result(6);
nmOff = result(7);
Lepi = result (8);

%%%%%%%%
% plot the data with fits
%%%

%%%%%
% model
%%
nm = nm0 + nmOff;
muaoxy = interp1(nmTG, muaoxyTG, nm);
muadeoxy = interp1(nmTG, muadeoxyTG, nm);
muawater = interp1(nmHQ, muaHQ, nm);
muamel = 6.6e11*nm.^-3.33;

%W = 0.65; % water content
%b = 1.0;  % scattering coefficient, musp = a*nm.^-b

Mie = 4.59e3*nm.^-0.913;
Ray = 1.74e12*nm.^-4;
musp = a*(Mie + Ray);

r = 0.2;
n = 1.4;

R = const*M;  % R = (Mskin/Mstd)(Rstd*Gstd/Gskin)

% epidermal melanin
Lepi = 2*0.0060;
Tepi = exp(-Mel*muamel*Lepi);

% +scatter only
mua = 1e-6;
pRS = const*getRdFarrell(mua, musp, n); 

% +scatter +Water
mua = W*muawater;
pRw  = const*getRdFarrell(mua, musp, n); 

% +scatter +Water +Mel
pRwM = const*Tepi.*getRdFarrell(mua, musp, n); 

% +scatter +Water +Mel +Blood
mua  = B*(S*muaoxy + (1-S)*muadeoxy) + W*muawater;
pRwMB = const*Tepi.*getRdFarrell(mua, musp, n); 
%%
% end model
%%%%%


figure(3);clf
set(figure(3),'position',[645    40   577   757],'color','w')
sz = 18;
subplot(2,1,1); hold off
plot(nm(u), M(u), 'ko','linewidth',1)
hold on
%plot(nm, pRS, 'k--','linewidth',2)
% plot(nm, pRw, 'b-','linewidth',2)
% plot(nm, pRwM, 'k-','linewidth',2)
plot(nm, pRwMB, 'r-','linewidth',2)
%
set(gca,'fontsize',sz,'linewidth',2)
xlabel('wavelength [nm]')
ylabel('M = M_s_k_i_n / M_s_t_d')
title('M = const T_e_p_i getRd(\mu_a, \mu_s'', n)')
axis([450 1700 0 1])
x = 410; ymax = 1; dy = .07;
text(x, ymax - dy, sprintf('W = %0.01f', W*100),'fontsize',sz)
text(x, ymax - 2*dy, sprintf('B = %0.4f', B),'fontsize',sz)
text(x, ymax - 3*dy, sprintf('S = %0.3f', S),'fontsize',sz)
text(x, ymax - 4*dy, sprintf('a = %0.3f', a),'fontsize',sz)
text(x, ymax - 5*dy, sprintf('const = %0.2f', const),'fontsize',sz)
text(x, ymax - 6*dy, sprintf('Mel = %0.3f', Mel),'fontsize',sz)
text(x, ymax - 7*dy, sprintf('Lepi = %0.3f',Lepi),'fontsize',sz)

subplot(2,1,2); hold off
semilogy(nm, W*muawater, 'b--','linewidth',2)
hold on
plot(nm, B*S*muaoxy, 'r-','linewidth',2)
plot(nm, B*(1-S)*muadeoxy, 'b-','linewidth',2)
plot(nm, Mel*muamel, 'k-','linewidth',2)
plot(nm, B*(S*muaoxy+(1-S)*muadeoxy)+W*muawater,'m-','linewidth',2)
plot(nm, Tepi,'k-','linewidth',2)
text(630, 1.4e-2,'total','fontsize',sz,'color','m')
text(480,0.07,'oxy','fontsize',sz,'color','r')
text(480,0.015,'deoxy','fontsize',sz,'color','b')
text(750,0.01,'water','fontsize',sz,'color','b')
text(650,0.5,'T.epidermal.melanin','fontsize',sz,'color','k')
set(gca,'fontsize',sz,'linewidth',2)
xlabel('wavelength [nm]')
ylabel('\mu_a [cm^-^1]')
axis([450 1700 1e-3 1])

figure(4);clf
set(figure(4),'position',[300    10   577   757],'color','w')
subplot(2,1,1); hold off
plot(nm, musp,'r-','linewidth',2)
hold on
plot(nm, a*Mie,'b-','linewidth',2)
plot(nm, a*Ray,'k-','linewidth',2)
text(600,140,['\mu_s''' sprintf(' = %0.2f(Mie + Ray)', a)],'fontsize',sz,'color','r')
text(600,120,'Mie = 4.x10^5 nm^-^1^.^5','fontsize',sz,'color','b')
text(600,100,'Ray = 2x10^1^2 nm^-^4','fontsize',sz,'color','k')
set(gca,'fontsize',sz)
xlabel('wavelength [nm]')
ylabel('\mu_s'' [cm^-^1]')
axis([450 1700 0 150])
%
subplot(2,1,2); hold off
%
nm2 = (100:1e4)';
Mie2 = 4.59e3*nm2.^-0.913;
Ray2 = 1.74e12*nm2.^-4;
musp2 = a*(Mie2 + Ray2);
%
loglog(nm, musp,'r-','linewidth',4,'color','r')
hold on
plot(nm, a*Mie,'b-','linewidth',4,'color','b')
plot(nm, a*Ray,'k-','linewidth',4,'color','k')
%
loglog(nm2, musp2,'r--','linewidth',2,'color','r')
plot(nm2, a*Mie2,'b--','linewidth',2,'color','b')
plot(nm2, a*Ray2,'k--','linewidth',2,'color','k')
%
text(600,400,['\mu_s''' sprintf(' = %0.2f(Mie + Ray)', a)],'fontsize',sz,'color','r')
%text(600,200,'Mie = 4.x10^5 nm^-^1^.^5','fontsize',sz,'color','b')
%text(600,100,'Ray = 2x10^1^2 nm^-^4','fontsize',sz,'color','k')
set(gca,'fontsize',sz)
xlabel('wavelength [nm]')
ylabel('\mu_s'' [cm^-^1]')
axis([100 1e4 1 1e4])

