function [h_recon, info] = recover_sparse_sequence_CP(f1, f2, start_point, param)

% Attempts to recover a sparse sequence whose selected DFT values are within limits
% using the Chambolle-Pock algorithm

primal = start_point;
% primal = start_point.primal;
dual = f2.L(start_point);

%step-sizes
% param.theta = 1;
% param.sigma = 1 / sqrt(f2.norm_L);
% param.tau = 1 / sqrt(f2.norm_L);

if param.debug_mode
    info.obj_eval = NaN(param.maxit,1);  %allocating for objective function history
end

info.time_algorithm = tic;

%% Algorithm
iter = 1; %iteration counter
x_n = primal;

while iter <= param.maxit
    if param.verbose > 1
        fprintf('Iteration no. %3i',iter)
%         fprintf('  relative change in coefficients:  %e \n', relat_change_coefs);
%         fprintf('  relative change in objective fun: %e \n', relat_change_obj);
        fprintf('\n');
    end
    
    %Using Moreau identity
%     dual = prox_adjoint( dual + param.sigma * f2.L(primal), param.sigma, f2);
    aux = dual + param.sigma * f2.L(primal);
    dual = aux - param.sigma * f2.prox(aux/param.sigma, 1/param.sigma);
    
    x_n_old = x_n;
    x_n = f1.prox( x_n - param.tau * f2.Lt(dual), param.tau);
    primal = x_n + param.theta * (x_n - x_n_old);
    
    if param.debug_mode
        info.obj_eval(iter) = f1.eval(primal) + f2.eval(f2.L(primal)); %record values of objective function
    end
    
    iter = iter + 1;
end

toc(info.time_algorithm)

sol.primal = primal;
sol.dual = dual;

info.iter = iter - 1;
disp(['Finished after ' num2str(info.iter) ' iterations.'])

h_recon = sol.primal; %output