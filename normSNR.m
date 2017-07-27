function [SNR, out] = normSNR(signal,noise,h,delay,L)
% Normalized SNR computed in the output of Blocking Matrix

[m, N] = size(signal);

if m>N
    signal = signal';
    noise = noise';
    [m, ~] = size(signal);
end
rnoise = residual(noise',h,delay);
rsignal = residual(signal',h,delay);

SNR = zeros(m,1);

out = zeros(size(signal));

for i = 1:m
    g = TDRTF(L,rnoise,noise(i,:)',delay,0.01);
    esig = filter(g,1,rsignal);
    enoi = filter(g,1,rnoise);
    out(i,:) = esig+enoi;
    SNR(i) = mean(esig.^2)/mean(enoi.^2);
end

SNR = mean(SNR);
