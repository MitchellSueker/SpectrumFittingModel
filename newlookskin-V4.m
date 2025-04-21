% lookskin1.m
% loads spectrum --> nm, skin, std, R
%
% M = (S T G D)/(S Tstd Gstd D) 
%   = T/Tstd*G/Gstd = T * constant
Rstd = 0.99; %REFLECTANCE PERCENTAGE

figure(1); clf
set(figure(1),'position',[31    42   580   755],'color','w')
sz = 18;

dinm = 3;

%%%%% Fig 1.1
subplot(2,1,1)

% Mskin = (table2array(readtable('Ref1.xlsx')));
% 
% nm = Mskin(:,1);
% skin = Mskin(:,2);
% %i=1; while nm(i)<400; i=i+1; end; i400 = i;
% %while nm(i)<900; i=i+1; end; i900 = i;
% %base = mean(skin(10:100));
% %skin = skin - base;
% plot(nm, skin,'r-','linewidth',2)
% hold on
% %
% 
% Mstd = (table2array(readtable('white.xlsx')));
% 
% std = Mstd(:,2);%*3;
% %base = mean(std(10:100));
% %std = std2- base;
% plot(nm, std,'b-','linewidth',2)
% set(gca,'fontsize',sz,'linewidth',2)
% xlabel('wavelength [nm]')
% ylabel('M')
% text(600,650,'Measurements','fontsize',sz)
% %axis([400 900 -100 700])
% legend('skin-dark','white -dark')


%%%%% Fig 1.2
% subplot(2,1,2)
data = (readtable('pantea2.xlsx'));
R = data{:, 2:end} %creates subset array
 nm = data{:, 1}
%% Apply First-Order Polynomial Smoothing (Savitzky-Golay Filter)
window_size = 11;  % Must be odd
poly_order = 1;  % First-order polynomial smoothing
R = sgolayfilt(R, poly_order, window_size);
% plot(nm, R,'r-','linewidth',2)
% set(gca,'fontsize',sz,'linewidth',2)
% xlabel('wavelength [nm]')
% ylabel('(Skin-dark)/(white -dark)')
% %('T = M T_s_t_d/GG')
% text(600, 0.9,'Reflectance','fontsize',sz)
%axis([400 900 0 1])

