function [h_recon, info] = recover_sparse_sequence_DR(f1, f2, start_point, param)
%
% Attempts to recover a sparse sequence whose selected DFT values are within limits
% using the Douglas-Rachford algorithm



%% Setting functions and params for Douglas-Rachford algorithm

DR_y = start_point;

info.time_algorithm = tic;

DR_x_old = DR_y;
         
relat_change_coefs = 1;
relat_change_obj = 1;
iter = 1; %iteration counter

if param.debug_mode
    info.obj_eval = NaN(param.maxit,1);  %allocating for objective function history
end


% while relat_change_coefs > 0.00001
% while relat_change_obj > param.tol
while iter <= param.maxit
    % DR: gamma = 1
    
    aux = f2.L(DR_y);
%     DR_x = f2.prox(DR_y,[]);
    DR_x = DR_y + f2.Lt( f2.prox(aux) - aux );
    
    if param.debug_mode
        info.obj_eval(iter) = f1.eval(DR_x) + f2.eval(DR_x); %record values of objective function
    end
    
    DR_y = DR_y + param.lambda*(f1.prox(2*DR_x-DR_y, param.gamma)-DR_x);
    if iter > 1
        if param.verbose > 1
            relat_change_coefs = norm(DR_x-DR_x_old) / norm(DR_x_old);
            relat_change_obj = norm(obj_eval(end) - obj_eval(end-1)) / norm(obj_eval(end-1));
            fprintf('Iteration no. %3i',iter)
            fprintf('  relative change in coefficients:  %e \n', relat_change_coefs);
            fprintf('  relative change in objective fun: %e \n', relat_change_obj);
            fprintf('\n');
        end
    end
    DR_x_old = DR_x;
    iter = iter + 1;
    
end

h_recon = DR_y;

toc(info.time_algorithm)

info.iter = iter - 1;
disp(['Finished after ' num2str(info.iter) ' iterations.'])

% info.obj_eval = obj_eval;