% function M = makeFig1Spectrum(name)
% The basic program for looking at one spectrum
% called [name '.Master.Scope']

function  A = getSpectrum(name)

file = strcat(name, '.Master.Scope');
%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(file ,'r');
for i=1:19; s = fgetl(fid); disp(s); end
A = fscanf(fid, '%f %f', [2 inf]);
A = A';
nm = A(:,1);
M = A(:,2);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%
