function [h_recon_CP, info] = CPsolver(H,S,mu)

N = length(H);

H = H / sqrt(N);  %same length, complex conjugated pairs present

subsample = find(S);   %randperm does not include 0, and it should not include the DC component (1)
sFh = H(S);

epsilons = mu(S); epsilons = epsilons(:);
    
f1.eval = @(x) norm(x,1);
f1.prox = @(x,T) sign(x).*max(abs(x)-T, 0);  %soft thresholding - works for complex numbers as well
f2.eval = @(x) eps; %simplified to zero (eps), since after each projection we are back at the constraint
f2.prox = @(x,T) epsilons .* (x-sFh) ./ max(abs(x-sFh),epsilons) + sFh; %complex box projection in the Fourier space
f2.L = @(x) subsampled_fft(x,subsample); %linear operator (complex); FFT is normalized to be unitary
f2.Lt = @(x) subsampled_fft_adjoint_and_conjugate(x,subsample,N); %its adjoint operator
f2.norm_L = 1; %spectral norm of the subsampled Fourier transform
start_point = zeros(N,1);
paramsolver.verbose = 1;  % display parameter
paramsolver.maxit = 600;        % maximum number of iterations
paramsolver.theta = 1;
paramsolver.sigma = 1 / sqrt(f2.norm_L);
paramsolver.tau = 1 / (paramsolver.sigma * f2.norm_L);
paramsolver.debug_mode = 1;  %debug mode
    


%%%%%%%% RUN THE ALGORITHM
[h_recon_CP, info] = recover_sparse_sequence_CP(f1, f2, start_point, paramsolver);
h_recon_CP = project_into_spectral_constraints(h_recon_CP, f2, subsample, N);
    

