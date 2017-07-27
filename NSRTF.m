function [h, H, VAR, Sxx, Syx, Svx] = NSRTF(x,y,NFFT,delay)
%

y = [zeros(delay,1); y(1:end-delay)];

win = hamming(NFFT);

X = stft(x',NFFT,NFFT/16,NFFT,win);
Y = stft(y',NFFT,NFFT/16,NFFT,win);

Sxx = transpose(X.*conj(X));
Syx = transpose(Y.*conj(X));


teta = zeros(2,NFFT/2+1);

for w = 1:NFFT/2+1 %225 Nfft/2+1 % 6:120% 16:64 % looping for each frequancy

  [teta(:,w)] = [Sxx(:,w) , ones(size(Sxx(:,w)))] \ Syx(:,w);

end

H   = [teta(1,1:NFFT/2+1) conj(teta(1,NFFT/2:-1:2))].'; % teta(1,:);%
Svx = [teta(2,1:NFFT/2+1) conj(teta(2,NFFT/2:-1:2))]; % teta(2,:);%

h = ifft(teta(1,:),NFFT,'symmetric')';

v = filter(h,1,x)-y;
V = stft(v',NFFT,NFFT/16,NFFT,win);
Svv = mean(V.*conj(V),2)';

VAR = real(Svv.*mean(1./Sxx)./(mean(Sxx).*mean(1./Sxx)-1))/length(x);
VAR = abs(VAR)';