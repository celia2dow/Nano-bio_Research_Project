function [l2_mean, l2_array] = l2_from_diffsMethod(times,av_data,tol,...
    freePrtcls_start_of_t, interactPrtcls_start_of_t,PARAMETERS,l1_mean)
% L2_FROM_DIFFSMETHOD estimates lambda2 by assuming that the gradient of
% the data (av_data) of interacting particles in a timestep before the CC 
% kicks in (times) is a constant when a threshold in the second derivative 
% is met (tol). Assume also that the gradient of associated particles in a 
% timestep will be equal to the number of free particles at the beginning 
% of that timestep (freePrtcls_start_of_t) multiplied by the rate of 
% transition lambda1 estimated via the method of differences (l1_mean). 
% Assume also that the gradient of internalised particles in a timestep 
% will be equal to the number of interacting particles at the beginning of 
% that timestep (interactPrtcls_start_of_t) multiplied by the dynamic lambda2.
%
% It returns an estimate per timestep as well as an overall mean
% calculated from the data beyond 0.2 hours for which the gradient of the
% data is approximately constant.
%
% Units of rates are per hour.

% Find when the steady state in interacting particles is reached, if it
% exists
deriv1_interacting = gradient(av_data(1,:),PARAMETERS.tstep_duration);
zero_deriv1_interacting_times = times(abs(deriv1_interacting(1:length(times)))<=tol);
zero_deriv1_interacting_times = zero_deriv1_interacting_times...
    (zero_deriv1_interacting_times > 0.2);
% Find the timesteps when a constant non-zero gradient in interacting
% particles is reached, if it exists
deriv2_interacting = gradient(deriv1_interacting,PARAMETERS.tstep_duration);
zero_deriv2_interacting_times = times(abs(deriv2_interacting(1:length(times)))<=tol);
zero_deriv2_interacting_times = zero_deriv2_interacting_times...
    (zero_deriv2_interacting_times > 0.2);

% If there is a turning point and thus a steady state - though note that
% this is never really appropriate
if any(zero_deriv1_interacting_times) && deriv1_interacting(end)<0
    % lambda_2 can be estimated by equating lambda_1*(# free particles) and
    % lambda_2*(# interacting particles) when a steady state is reached in
    % the number of iteracting particles per cell
    indices = zeros(1,length(zero_deriv1_interacting_times));
    for i = 1:length(zero_deriv1_interacting_times)
        indices(i) = find(times == zero_deriv1_interacting_times(i));
    end
    indices = indices -1;
    % Estimate lambda_2 using 
    %       lambda_1*(# free particles at start of t) = 
    %       lambda_2* (# interacting particles at start of t) 
    l2_array = l1_mean .* freePrtcls_start_of_t(indices) ./ ...
        interactPrtcls_start_of_t(indices);
    l2_mean = mean(l2_array);
    
% If there are points of constant gradient
elseif zero_deriv2_interacting_times 
    % lambda_2 can be estimated by equating the gradient of the curve of 
    % average interacting particles over time at a timestep with the 
    % difference in gradients of the curve for average associated particles 
    % over time and that of internalised particles over time. This is done
    % at times when the interacting gradient is found to be approximately
    % constant.
    indices = zeros(1,length(zero_deriv2_interacting_times));
    for i = 1:length(zero_deriv2_interacting_times)
        indices(i) = find(times == zero_deriv2_interacting_times(i));
    end
    indices = indices -1;
    % Find the average gradient at each timestep. Estimate lambda_2 using 
    %       (gradient of interacting curve) * tstep_duration = 
    %       lambda_1*(# free particles at start of t) - 
    %       lambda_2* (# interacting particles at start of t)
    l2_array = (l1_mean .* freePrtcls_start_of_t ...
            - deriv1_interacting(2:end))./...
            interactPrtcls_start_of_t;
    l2_mean = mean(l2_array(indices));
else
    fprintf("\nlambda2 was not able to be estimated via differences due to noise \n")
    l2_array = [];
    l2_mean = [];
end