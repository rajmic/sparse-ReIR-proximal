function experimentSOCP(intlength)
% intlength ... the number of samples used within each trial

interval = 1:intlength; % processed batch of data
Fs = 16000;
FFTSHIFT = 256;
L = 1024; % FFT length 
D = 20; % global delay

ntrials = 100; % # trials

% Directory of Multichannel Impulse Response Database
% Available at http://www.eng.biu.ac.il/gannot/downloads/
pathData = 'd:\Databaze\Impulse\';
% parameters for the database
micSpac = '3-3-3-8-3-3-3'; % microphone-array spacing
dist = 2; % source distance [m]
T60 = 610; %160 360 610
micsUsed = [3 4]; % used microphones

iniSNR = 0; % initial SNR

% Path to a WAV (speech) file
soundPath = 'dev_Sq1_Co_B_src.wav';

profile on

for trial = 1:ntrials

% Speaker Spatial Images    
angle = 0;
a = SP.loadImpulseResponse(pathData, micSpac, angle, dist, T60, 0.6);
a.resample_it(Fs);a = a.to_double();a = a(:,micsUsed);
s1 = SP.signal(soundPath);s1.resample_it(Fs);s1 = s1.to_double();s1=s1/norm(s1)*10;
ofset = floor((length(s1)-intlength)*rand); 
s1 = s1(ofset:ofset+intlength-1); % random block of length given by intlength
resp = zeros(2,length(s1),2);
for lp = 1:2, resp(lp,:,1) = filter(a(:,lp),1,s1)'; end

% Noise Spatial Images
angle = 315;
a = SP.loadImpulseResponse(pathData, micSpac, angle, dist, T60, 0.6);
a.resample_it(Fs);a = a.to_double();a = a(:,micsUsed);
s1 = randn(size(s1)); for lp=1:2, resp(lp,:,2) = filter(a(:,lp),1,s1)';end
% Setting initial SNR
resp(:,:,2)=resp(:,:,2)*(sqrt((mean(mean(resp(:,:,1).^2))/mean(mean(resp(:,:,2).^2))))*10^(-iniSNR/20));

% mixed signals
x = sum(resp,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison of various RTF estimates

%%% Noise-free estimates

% Time domain noise-free RTF (ReIR) estimate
h = TDRTF(L,resp(1,interval,1)',resp(2,interval,1)',D);
oSNR_T=10*log10(normSNR(resp(:,interval,1),resp(:,interval,2),h,D,L)); 

% Frequency domain noise-free RTF estimate
g = FDRTF(L,resp(1,interval,1)',resp(2,interval,1)',D,32);
oSNR_F=10*log10(normSNR(resp(:,interval,1),resp(:,interval,2),g,D,L)); 

% Frequency domain nonstationarity-based noise-free RTF estimate
h2 = NSRTF(resp(1,interval,1)',resp(2,interval,1)',L,D);
oSNR_F2 = 10*log10(normSNR(resp(:,interval,1),resp(:,interval,2),h2,D,L));

%%% Noisy RTF estimates

% Frequency domain RTF estimate
g2 = FDRTF(L,x(1,interval)',x(2,interval)',D,32);
g2 = g2(1:L);
preSNR = 10*log10(normSNR(resp(:,interval,1),resp(:,interval,2),g2,D,L));

% Nonstationarity-based noisy RTF estimate
[g4, G4, nsrtf_var] = NSRTF(x(1,interval)',x(2,interval)',L,D);
preSNR3 = 10*log10(normSNR(resp(:,interval,1),resp(:,interval,2),g4,D,L));%10*log10(mean(residual(resp(:,interval,1)',g2,D).^2)/mean(residual(resp(:,interval,2)',g2,D).^2));


% Selection of ``reliable'' frequencies using the first channel as a
% reference

XL = stft(x(1,interval),L,FFTSHIFT);

% kurtosis
kurt = (mean(abs(XL).^4,2)-abs(mean(XL.^2,2)).^2)./mean(abs(XL).^2,2).^2 - 2; 

% frequency-dependent SNR
A = stft(resp(1,interval,1),L,FFTSHIFT,L,ones(L,1));
B = stft(resp(1,interval,2),L,FFTSHIFT,L,ones(L,1));
localSNR = sum(abs(A).^2,2)./sum(abs(B).^2,2);
clear A B

% oracle selection
[~, orderoracle] = sort(localSNR(1:L/2+1),'descend');
% kurtosis-based selection
[~, orderkurtosis] = sort(kurt(1:L/2+1),'descend');

sparirtime = 0;
CPtime = 0;
DRtime = 0;

percs = [2:2:20 25:5:100];
index_p = 0;
for percentage = percs
    index_p = index_p+1;

    fprintf('Percentage: %d, trial: %d\n',percentage,trial)

    % Oracle selection
    S = false(L,1);
    S(orderoracle(1:floor((L/2+1)/100*percentage))) = true;

    % SpaRIR (LASSO)
    tic
    [a, NumIt] = SpaRIR(G4,S,D);
    sparirtime = sparirtime + toc/NumIt;
    SNR_oracle_WLASSO(trial,index_p) = normSNR(resp(:,interval,1),resp(:,interval,2),a,D,L); 
    
    % CP (SOCP)
    tic
    [a, info] = CPsolver(G4, S, 2*sqrt(nsrtf_var));
    CPtime = CPtime + toc/600;
    SNR_oracle_CP(trial,index_p) = normSNR(resp(:,interval,1),resp(:,interval,2),a,D,L); 

    % DR (SOCP)
    tic
    a = DRsolver(G4, S, 2*sqrt(nsrtf_var));
    DRtime = DRtime + toc/600;
    SNR_oracle_DR(trial,index_p) = normSNR(resp(:,interval,1),resp(:,interval,2),a,D,L);
    
    % Kurtosis-based selection
    
    S = false(L,1);
    S(orderkurtosis(1:floor((L/2+1)/100*percentage))) = true;

    % SpaRIR (LASSO)
    [a, NumIt] = SpaRIR(G4,S,D);
    SNR_kurt_WLASSO(trial,index_p) = normSNR(resp(:,interval,1),resp(:,interval,2),a,D,L);
    
    % CP (SOCP) 
    [a, info] = CPsolver(G4, S, 2*sqrt(nsrtf_var));
    SNR_kurt_CP(trial,index_p) = normSNR(resp(:,interval,1),resp(:,interval,2),a,D,L);

    % DR (SOCP) - for comparison, controlled purely through the selection of
    % epsilon's (not directly through the set S)
    aux = ones(L,1);
    aux(S) = sqrt(nsrtf_var(S));
    SS = false(L,1);
    SS(1:L/2+1) = true;
    a = DRsolver(G4, SS, aux);
    SNR_kurt_DR(trial,index_p) = normSNR(resp(:,interval,1),resp(:,interval,2),a,D,L);
    
end
end
profile report

figure
plot(percs,10*log10(mean(SNR_oracle_WLASSO,1)))
hold on
plot(percs,10*log10(mean(SNR_oracle_CP,1)))
plot(percs,10*log10(mean(SNR_oracle_DR,1)))
plot(percs,10*log10(mean(SNR_kurt_WLASSO,1)),'--')
plot(percs,10*log10(mean(SNR_kurt_CP,1)),'--')
plot(percs,10*log10(mean(SNR_kurt_DR,1)),'--')
hold off
xlabel('percentage [%]')
ylabel('normSNR [dB]')
