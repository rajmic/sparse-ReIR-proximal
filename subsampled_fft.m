function y = subsampled_fft(x,indexes)
    % performs fft and subsamples it
    % the FFT is normalized (to become unitary)
    aux = fft(x) / sqrt(length(x));
    y = aux(indexes);
end