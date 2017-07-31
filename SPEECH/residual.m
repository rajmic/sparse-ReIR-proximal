function [res, FB]=residual(x,h,delay)

res=filter(h,1,x(:,1))-[zeros(delay,1);x(1:end-delay,2)];
FB=(filter(h,1,x(:,2))+[zeros(delay,1);x(1:end-delay,1)])/2;
%fprintf('Original RMS [dB]: %2.2f\n',10*log10(mean(x(:).^2)));
%fprintf('Residual variance [dB]: %2.2f\n',10*log10(mean(res.^2)));
