function out = project_into_spectral_constraints(in, f2, subsample, N);
%project the time signal into spectral constraints
projected = f2.prox( f2.L(in) );
aux = fft(in) / sqrt(N); %original spectrum
aux(subsample) = projected; %replace the important entries with the projections
aux(N-subsample+2) = conj(projected); %symmetric complex conjugates
aux = aux(1:N);  %in case of presence of DC, cut its reflection
aux(1) = real(aux(1));  %DC component must be real (projection)
aux(N/2+1) = real(aux(N/2+1));  %component in the middle must be real (projection)
out = ifft(aux) * sqrt(N);
end