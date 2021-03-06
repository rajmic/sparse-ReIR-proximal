function [g, iter] = SpaRIR(G,S,delay,weights,q,gini)
% Sparse approximation of (relative) impulse response from incomplete 
% measurements of the corresponding (relative) transfer function 
%
% USAGE: g = SpaRIR(G,S,delay,weights,q,gini)
%
% inputs:
%
% G ... (incomplete) transfer function; an Lx1 vector; L is the length of g
% S ... an Lx1 logical array such that G(S) is the incomplete measurement 
%       of the transfer function; the values of G(~S) are not used
% delay ... expected index of the main peak in g; default value = 20; this
%           value need not have any influence on the result
% weights ... an Lx1 weighting vector or scalar; the higher the value of 
%             weights(i), the higher the probability that g(i) is zero; 
%             if scalar, all weights are the same; if empty, default value
%             is used
% q ... oversampling factor (default q=1)
% gini ... initialization for the computation of g
%
% outputs:
%
% g ... the sparse approximation of the impulse response (real-valued)
%
% Coded by Zbynek Koldovsky, last change January 2016
%
% Please cite in case of use:
%
% Z. Koldovsky, J. Malek and S. Gannot, "Spatial Source Subtraction based 
% on Incomplete Measurements of Relative Transfer Function," submitted for
% publication in 2014.
%

% Author(s): Zbyn�k Koldovsk�
% Technical University of Liberec
% Studentsk� 1402/2, LIBEREC
% Czech Republic
%
%
% This is unpublished proprietary source code of TECHNICAL UNIVERSITY OF
% LIBEREC, CZECH REPUBLIC.
% 
% The purpose of this software is the dissemination of scientific work for
% scientific use. The commercial distribution or use of this source code is
% prohibited. The copyright notice does not evidence any actual or intended
% publication of this code. Term and termination:
% 
% This license shall continue for as long as you use the software. However,
% it will terminate if you fail to comply with any of its terms and
% conditions. You agree, upon termination, to discontinue using, and to
% destroy all copies of, the software.� Redistribution and use in source and
% binary forms, with or without modification, are permitted provided that
% the following conditions are met:
% 
% Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer (Limitation of
% warranties and liability). Redistributions in binary form must reproduce
% the above copyright notice, this list of conditions and the following
% disclaimer (Limitation of warranties and liability) in the documentation
% and/or other materials provided with the distribution. Neither name of
% copyright holders nor the names of its contributors may be used to endorse
% or promote products derived from this software without specific prior
% written permission.
% 
% The limitations of warranties and liability set out below shall continue
% in force even after any termination.
% 
% Limitation of warranties and liability:
% 
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE HEREBY
% DISCLAIMED.� IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS� OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
% OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
% SUCH DAMAGE.

L = length(G);

if rem(L,2)==1, error('The length of G must be even.'); end  

S=S(1:L/2+1); % S(k) and G(k) for k > L/2+1 not used due to the symmetry of FFT (g is assumed to be real-valued)

if ~exist('q','var')
    q=1;
end

if q>1
    G=[G(1:L/2+1);zeros(q*L-L/2-1,1)];
    S=[S;false(q*L-L/2-1,1)];
    L = q*L;
end

        
y=[real(G(S));imag(G(S))];
M=length(y);

if ~exist('delay','var')
    delay = 20;
end

delay = q*delay;

if ~exist('gini','var') % if no initialization is given
    g = zeros(L,1);    
    g(delay+1)=1;
else
    g = gini;
end


if ~exist('weights','var')
    tau = sqrt(L)/(y'*y);% * ones(L,1);
    tau = tau*exp(0.11*abs(((1:L)'-delay)/q).^0.3);
elseif isempty(weights)
    tau = sqrt(L)/(y'*y);
    tau = tau*exp(0.11*abs(((1:L)'-delay)/q).^0.3);
elseif length(weights)==1
    tau = ones(L,1)*weights;
else
    tau = reshape(repmat(weights',q,1),L,1);
end

maxiter = 10000;
alphamax = 1e5; % maximum step-length parameter alpha
alphamin = 1e-7; % minimum step-length parameter alpha
tol=1e-4;

aux=zeros(L/2,1);
G = fft(g);
Ag = [real(G(S));imag(G(S))];
r = Ag - y; % instead of r = A*g - y;
aux(S(1:L/2+1))=r(1:M/2) + 1i*r(M/2+1:end);
gradq=L/2*ifft(aux,L,'symmetric'); % instead of gradq = A'*r;
alpha = 10;
support = g~=0;
iter = 0;

crit = zeros(maxiter,1);
criterion = -tau(support).*sign(g(support))-gradq(support);
crit(iter+1) = sum(criterion.^2);


while (crit(iter+1) > tol) && iter<maxiter
    prev_r = r;
    prev_g = g;
    g = soft(prev_g - gradq*(1/alpha),tau/alpha);
    dg = g - prev_g;
    DG = fft(dg);
    Adg = [real(DG(S));imag(DG(S))];
    r = prev_r + Adg; %faster than A*g - y;
   	dd  = dg(:)'*dg(:);  
	dGd = Adg(:)'*Adg(:);
	alpha = min(alphamax,max(alphamin,dGd/(realmin+dd)));
    iter = iter + 1;
    support = g~=0;
    aux(S(1:L/2+1))=r(1:M/2) + 1i*r(M/2+1:end);
    gradq=L/2*ifft(aux,L,'symmetric');
    criterion = -tau(support).*sign(g(support))-gradq(support);
    crit(iter+1) = sum(criterion.^2) + sum(abs(gradq(~support))-tau(~support)>tol);
end

if q>1
    g=q*resample(g,1,q,100);
end

fprintf('SpaRIR: %d iterations done.\n',iter);

if sum(g==0)==L
    disp('\n Zero solution! Do you really like it?') 
    disp('Check input parameters such as weights (tau), delay, gini or')
    disp('internal parameters such as tol, alphamax, alphamin, etc.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y = soft(x,T)
if sum(abs(T(:)))==0
   y = x;
else
   y = max(abs(x) - T, 0);
   %y = sign(x) .* y;
   y = y./(y+T) .* x;
end