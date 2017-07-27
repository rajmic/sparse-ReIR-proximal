function X = stft(x, n, wshift, NFFT, window)

if(nargin<4)
   NFFT = n;  
end    

if(nargin<5)
   window = boxcar(n);  
end    

N = length(x);

N = ceil((N-n)/wshift)*wshift + n;

x = [x zeros(1,N-length(x))];

nwindow=(N-n)/wshift;

X=zeros(NFFT, nwindow + 1);

for i=0:nwindow
    X(1:n,i+1)=x(i*wshift+1:i*wshift+n)'.*window;
end
X=fft(X);