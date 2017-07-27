function x = istft(X, wshift, NFFT, window)

[n m] = size(X);

if nargin<3
    NFFT=n;
end

if(nargin<4)
   window = boxcar(n);  
end  

Y=real(ifft(X,NFFT,'symmetric'));
x=zeros(1,(m-1)*wshift+NFFT);
 
for i=0:m-1
    x(i*wshift+1:i*wshift+NFFT)=x(i*wshift+1:i*wshift+NFFT)+Y(1:NFFT,i+1)'.*window';
end

k=NFFT/wshift;
for i=2:k
    x((i-1)*wshift+1:i*wshift)=x((i-1)*wshift+1:i*wshift)/i;
end

x(k*wshift+1:m*wshift)=x(k*wshift+1:m*wshift)/k;

for i=m:m+k-3
    x(i*wshift+1:(i+1)*wshift)=x(i*wshift+1:(i+1)*wshift)/(m+k-1-i);
end
