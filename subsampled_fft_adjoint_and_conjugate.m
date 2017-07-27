function y = subsampled_fft_adjoint_and_conjugate(x,indexes,N)
    % Conjugates the complex pairs (can be seen as a projection) and performs adjoint operation to 'subsampled_fft'
    % Lower "half" of the spectrum is supposed to be established
    %
    % The computation is done by taking IFFT of half spectrum and taking real part.
    % The IFFT is turned to unitary by normalization.

    aux = zeros(N,1); %prepare array
    aux(indexes) = 2*x; %input values
    
    %the following will have an impact only if the freequencies are
    %indexed; otherwise it works with zeros
    aux(1) = real(aux(1)/2);  %DC component must be real (projection)
    aux(N/2+1) = real(aux(N/2+1)/2);  %component in the middle must be real (projection)

    y = real(ifft(aux) * sqrt(N));

%     aaux = zeros(N,1); %prepare array
%     aaux(indexes) = x; %input values
%     aaux(N-indexes+2) = conj(x); %symmetric complex conjugates
%     aaux = aaux(1:N); %if DC component is present, cut its copy
%     
%     aaux(1) = real(aaux(1));  %DC component must be real (projection)
%     aaux(N/2+1) = real(aaux(N/2+1));  %component in the middle must be real (projection)
%     
%     ay = ifft(aaux) * sqrt(N); %should be real
%     
%     if norm(y-ay) > 10e-13
%         error('norms!')
%     end
end