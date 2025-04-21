function Rd = getRdFarrell(mua, musp, n)
% function Rd = getRdFarrell(mua, musp, n)
%  	n = ratio n_tissue/n_air
%	mua = absorption coefficient, eg., [cm^-1] or [mm^-1]
%	musp = reduced scattering coefficient eg., [cm^-1] or [mm^-1]
%   EITHER mua and musp are vectors, and r is scaler
%   OR r is vector and mua,musp are scalars
%   n is always a scaler = ntissue/nmedium.outside.
% 	Farrell TJ, MS Patterson, B Wilson,
%  	Medical Physics 19:879-888, 1992
% slj 9/2009

ri = 0.6681 + 0.0636*n + 0.7099./n - 1.4399./n.^2;
A = (1 + ri)./(1 - ri);

zo = 1./(mua + musp);
D = zo/3;
delta = sqrt(D./mua);
mueff = 1./delta;
ap = musp./(mua+ musp);

Rd = ap.*exp(-mueff.*zo)/2.*(1 + exp(-4/3*A.*sqrt(3*(1-ap))));




