function [BPM, s] = SSR(y, Fs, N)
if (nargin == 2)
% resolution is set to 1 BPM
N = 60*Fs;
end
M = length(y);
% construct fourier matrix
Phi = zeros(M,N);
complex_factor = 1i*2*pi/N;
for m = 1:M
for n = 1:N
Phi(m,n) = exp(complex_factor*(m-1)*(n-1));
end
end
% calculat
BPM = 60*Fs/N*((1:N)-1);
% construct sparse spectrum
s = FOCUSS(y,Phi);
% Band-pass frequencies in BPM
BPM_bp = [40 200];
% Band-pass spectrum
% s = BP(s, Fs, BPM_bp);
 % function [ SSRsig ] = BP( SSRsig, Fs, BPM_bp )

N = length(s);
f_lo = BPM_bp(1) / 60;
f_hi = BPM_bp(2) / 60;
R = floor(f_lo/Fs*N+1) : ceil(f_hi/Fs*N+1);
H = zeros(N,1);
H(R) = 1;
s = s .* H;
end