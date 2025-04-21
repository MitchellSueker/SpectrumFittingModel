function err = fitspectrum2(start, nm0, M, nmTG, muaoxyTG, muadeoxyTG, nmHQ, muaHQ, u, v)
global cnt

cnt = cnt + 1;
GRAPH = 1; % 1 = ON, 0 = OFF
W = start(1); % water content

B     = start(2);
S     = start(3);
a     = start(4);
Mel   = start(5);
const = start(6);
nmOff = start(7);

%%%%%
% model
%%
nm = nm0 + nmOff;
muaoxy = interp1(nmTG, muaoxyTG, nm);
muadeoxy = interp1(nmTG, muadeoxyTG, nm);
muawater = interp1(nmHQ, muaHQ, nm);
muamel = 6.6e11*nm.^-3.33;
%figure; plot(muaoxy)

% Mitch: check the scattering model used for skin (unit in cm)

Mie = 4.59e3*nm.^-0.913;
Ray = 1.74e12*nm.^-4;
musp = a*(Mie + Ray);

r = 0.0;
n = 1.4;

% epidermal melanin
Lepi = 2*start(8);
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

err1 = sum( (pRwMB - M).^2 );
err2 = sum( ((pRwMB(v) - M(v))).^2);
err  = err1*err2; 
err = err*(max(pRS)-const > 0.02);

if S>1;   err = err*10; end
if S<0.0; err = err*10; end

if cnt<500
	FLAG = (mod(cnt,10)==0);
else
	FLAG = (mod(cnt,100)==0);
end

if ( GRAPH & FLAG )
figure(2); hold off
set(figure(2),'position',[645    40   577   757],'color','w')
sz = 18;
subplot(2,1,1); hold off
plot(nm, M, 'ko')
hold on
plot(nm(v), M(v), 'go')
plot(nm, pRS, 'k--','linewidth',2)
plot(nm, pRw, 'b-','linewidth',2)
plot(nm, pRwM, 'k-','linewidth',2)
plot(nm, pRwMB, 'r-','linewidth',2)

set(gca,'fontsize',sz)
xlabel('wavelength [nm]')
ylabel('M = M_s_k_i_n/M_s_t_d')
axis([450 1700 0 1])
x = 410; ymax = 1; dy = .07;
text(x, ymax - dy, sprintf('W = %0.01f', W*100),'fontsize',sz)
text(x, ymax - 2*dy, sprintf('B = %0.4f', B),'fontsize',sz)
text(x, ymax - 3*dy, sprintf('S = %0.3f', S),'fontsize',sz)
text(x, ymax - 4*dy, sprintf('a = %0.3f', a),'fontsize',sz)
text(x, ymax - 5*dy, sprintf('const = %0.2f', const),'fontsize',sz)
text(x, ymax - 6*dy, sprintf('Mel = %0.3f', Mel),'fontsize',sz)
text(x, ymax - 7*dy, sprintf('nmOff = %0.3f',nmOff),'fontsize',sz)
drawnow

end


