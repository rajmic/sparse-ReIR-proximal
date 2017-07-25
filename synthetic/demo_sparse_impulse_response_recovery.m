function demo_sparse_impulse_response_recovery

% written by Pavel Rajmic, Brno University of Technology, 2017

clc
clear
close all


%% Parameters
N = 128; %length of the signal; N must be even (because of complex conjugation)
k = 5; %sparsity of the signal in the time domain

%percentage of known spectral coefficients (within the half-spectrum, i.e. coefs 1,...,N/2)
percent_spect = 30; 

epsilon = .05; %tolerance in the modulus spectrum (will be randomized entrywise below)

%set colors etc. for plotting results
color_orig = 'b';
color_CP = 'r';
color_DR = 'g';
barWidth_orig = .95;
barWidth_CP = .65;
barWidth_DR = .25;

%% Generated (original) sparse signal
% dictionary in which the signal will be sparse
A = @(x) (x);  %Identity dictionary, i.e. signal is sparse itself

if k>N
    error('Sparsity cannot be greater than signal length');
end

% generate random sparse support
support = randperm(N)';  %randperm does not include 0
support = sort(support(1:k));

% Determine values of the coefficients
x_support =  1 + 3*rand([k 1]);
x_support =  x_support .* (((rand([k 1])>.5)*2)-1); %randomizing signs
x = zeros(N,1);
x(support) = x_support; %complete vector including zero coefs

% synthesize signal
h = A(x);  %target impulse response


%% Compute spectrum and subsample it
Fh = fft(h) / sqrt(N);  %same length, complex conjugated pairs present
% random subsampling
% subsample = (randperm(N/2+1))';   %randperm(x) permutes 1:x;
subsample = (randperm(N/2)+1)';   %randperm(x) permutes 1:x; do not include the DC component (1) !
ks = round(percent_spect*N/100/2);  %how many spectral samples will be taken
subsample = sort(subsample(1:ks)) %a fraction of FFT values, sorted (first half-spectrum)
% subsample = ...
% [     7
%     16
%     22
%     23
%     25
%     27
%     28
%     29
%     30
%     32
% ]

sFh = Fh(subsample);  %subsampled spectrum

%create spectral tolerances
epsilons = epsilon * ones(ks,1); %make a vector out of a single value
epsilons = epsilons + (0.03 * rand(ks,1) - 0.005 ); %randomize it 

%or just load a fixed dataset
% clear
% load('imp_response_dataset_04')

% show impulse response
fig_coef = figure;
handle_coefs_orig = bar(h, color_orig);
set(handle_coefs_orig,'BarWidth', barWidth_orig)
hold on;
axis tight
title('Original coefficients of sparse signal');

% show spectrum and limits
fig_spectrum = figure;
plot(subsample,epsilons,'b*')
hold on
title('Magnitude difference between (half-)spectra of original and reconstructed signals');
axis tight


%% Chambolle-Pock algorithm
disp('%%%%%%%%%%%%%%% Chambolle-Pock algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%l1-norm for sparsity in time domain
f1.eval = @(x) norm(x,1);
f1.prox = @(x,T) sign(x).*max(abs(x)-T, 0);  %soft thresholding - works for complex numbers as well

%indicator function of the feasible set
f2.eval = @(x) eps; %simplified to zero (eps), since after each projection we are back at the constraint
% f2.eval = @(x) 1 / (1 - ( any( abs(x(subsample)-sFh ) >= epsilon ))) - 1; %zero if x is in the feasible set, Inf otherwise
f2.prox = @(x,T) epsilons .* (x-sFh) ./ max(abs(x-sFh),epsilons) + sFh; %complex box projection in the Fourier space
% f2.prox = @(x,T) proj_box_complex(x-sFh,epsilons) + sFh; %complex box projection in the Fourier space
f2.L = @(x) subsampled_fft(x,subsample); %linear operator (complex); FFT is normalized to be unitary
f2.Lt = @(x) subsampled_fft_adjoint_and_conjugate(x,subsample,N); %its adjoint operator
% f2.Lt = @(x) subsampled_fft_adjoint(x,subsample,N); %its adjoint operator
f2.norm_L = 1; %spectral norm of the subsampled Fourier transform

% start_point = subsampled_fft_adjoint_and_conjugate(sFh,subsample,N);
% start_point = subsampled_fft_adjoint(sFh,subsample,N);
start_point = zeros(N,1);

% paramsolver.verbose = paramsolver.verbose;
% paramsolver.maxit = paramsolver.maxit;
% % paramsolver.nu = paramsolver.nu;
% paramsolver.theta = 1;
% paramsolver.sigma = 1 / sqrt(f2.norm_L);
% paramsolver.tau = 1 / sqrt(f2.norm_L);
% paramsolver.tau = paramsolver.tau;
% paramsolver.sigma = paramsolver.sigma;
% paramsolver.tol = paramsolver.tol;
% paramsolver.debug_mode = paramsolver.debug_mode;

% setting different parameter for the simulation
paramsolver.verbose = 1;  % display parameter
paramsolver.maxit = 150;        % maximum number of iterations
% paramsolver.nu = 1;       % spectral norm of subsampled FFT (FFT is normalized to be unitary)
paramsolver.theta = 1;
paramsolver.sigma = 1 / sqrt(f2.norm_L);
% paramsolver.tau = 1 / sqrt(f2.norm_L);
paramsolver.tau = 1 / (paramsolver.sigma * f2.norm_L);
% paramsolver.tau = 0.999;       % primal stepsize
% paramsolver.sigma = 0.999;       % dual stepsize
% paramsolver.rescale = true;      % use rescaled variant of algorithm (see Komodakis)
% paramsolver.tol = 1e-8;        % Tolerance to stop iterating
paramsolver.debug_mode = 1  %debug mode


%%%%%%%% RUN THE ALGORITHM
% time_algorithm = tic;
% [h_recon, info] = chambolle_pock(start_point, f1, f2, paramsolver);
% toc(time_algorithm)
[h_recon_CP, info_CP] = recover_sparse_sequence_CP(f1, f2, start_point, paramsolver);
% [h_recon, info, timing] = recover_sparse_sequence_CP(f1, f2, start_point, paramsolver);

% final projection into the spectral constraints
h_recon_CP = project_into_spectral_constraints(h_recon_CP, f2, subsample, N);


% Display results
% timing
complexness_ratio = norm(imag(h_recon_CP)) / norm(real(h_recon_CP));  %check realness
disp(['Ratio of norms, imaginary part vs. real part: ' num2str(complexness_ratio)]);

disp('Distance of prescribed and computed spectral coefficients (before projection), and the limits:')
H_recon_CP = fft(h_recon_CP)/sqrt(N);
abs_difference_CP = abs( H_recon_CP(subsample) - sFh );
[abs_difference_CP epsilons]

% h_recon = ifft( f2.prox(temp,[]) ) * sqrt(N);
% % disp(['Distance of prescribed and computed spectral coefficients: 'num2str(___________)]);
% show (difference) spectrum

% show spectrum
figure(fig_spectrum)
handle_spectrum_recon_CP = bar(subsample, abs_difference_CP, color_CP);
set(handle_spectrum_recon_CP, 'BarWidth', barWidth_CP)
%limits
% hold on
% plot(subsample,epsilons,'b*')
% title('Magnitude difference between (half-)spectra of original and reconstructed signals');
axis tight

% show impulse response
figure(fig_coef);
handle_coefs_recon_CP = bar(real(h_recon_CP), color_CP);
set(handle_coefs_recon_CP,'BarWidth',barWidth_CP)
title('Original and reconstructed coefficients (real part) of sparse signal');
axis tight

% show convergence plot
if paramsolver.debug_mode
    fig_obj = figure;
%     figure(fig_obj);
    handle_obj_CP = plot(info_CP.obj_eval, color_CP);
%     set(handle_coefs_recon_CP,'BarWidth',0.5)
    axis tight
    hold on
    title('Objective function w.r.t. iterations');
end


%% Douglas-Rachford (composite) algorithm
disp('%%%%%%%%%%%%%%% Douglas-Rachford %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

f1 = f1;
f2 = f2;
% clear f3
start_point = start_point;

% setting different parameter for the simulation
paramsolver.verbose = paramsolver.verbose;
paramsolver.lambda = 1;   % step for DR algorithm (default 1)
paramsolver.gamma = 1 ; %1e-1;   % in this problem gamma plays the role of threshold for soft thresholding
% paramsolver.tol = paramsolver.tol;

paramsolver

%%%%%%%% RUN THE ALGORITHM
[h_recon_DR, info_DR] = recover_sparse_sequence_DR(f1, f2, start_point, paramsolver);

% final projection into the spectral constraints
h_recon_DR = project_into_spectral_constraints(h_recon_DR, f2, subsample, N);

% Display results
complexness_ratio = norm(imag(h_recon_DR)) / norm(real(h_recon_DR));  %check realness
disp(['Ratio of norms, imaginary part vs. real part: ' num2str(complexness_ratio)]);

disp('Distance of prescribed and computed spectral coefficients (before projection), and the limits:')
H_recon_DR = f2.L(h_recon_DR);
abs_difference_DR = abs( H_recon_DR - sFh );
[abs_difference_DR epsilons]

% show convergence plot
%%% beware that the objective fun is evaluated before the projection into spectral constraints
if paramsolver.debug_mode
    figure(fig_obj);
    handle_obj_DR = plot(info_DR.obj_eval, color_DR);
%     set(handle_coefs_recon_CP,'BarWidth',0.5)
    axis tight
    legend([handle_obj_CP handle_obj_DR],'CP', 'DR')
end

% show spectrum
figure(fig_spectrum)
handle_spectrum_recon_DR = bar(subsample, abs_difference_DR, color_DR);
set(handle_spectrum_recon_DR,'BarWidth', barWidth_DR)
% %limits
% hold on
handle_epsilons = plot(subsample,epsilons,'b*');
% title('Magnitude difference between (half-)spectra of original and reconstructed signals');
axis tight
legend([handle_epsilons handle_spectrum_recon_CP handle_spectrum_recon_DR],'limits', 'CP', 'DR')

% show impulse response
figure(fig_coef);
handle_coefs_recon_DR = bar(real(h_recon_DR), color_DR);
set(handle_coefs_recon_DR,'BarWidth', barWidth_DR)
% title('Original and reconstructed coefficients (real part) of sparse signal');
axis tight
legend([handle_coefs_orig handle_coefs_recon_CP handle_coefs_recon_DR],'original', 'CP', 'DR')

end