function [l2_mean, l2_array] = l2_from_diffsMethod(times,av_data,tol,...
    freePrtcls_start_of_t, interactPrtcls_start_of_t,PARAMETERS,l1)
% L2_FROM_DIFFSMETHOD estimates lambda2 by assuming that the gradient of
% the data (av_data) of interacting particles in a timestep before the CC 
% kicks in (times) is a constant when a threshold in the second derivative 
% is met (tol). Assume also that the gradient of associated particles in a 
% timestep will be equal to the number of free particles at the beginning 
% of that timestep (freePrtcls_start_of_t) multiplied by the rate of 
% transition lambda1 estimated via the method of differences (l1). 
% Assume also that the gradient of internalised particles in a timestep 
% will be equal to the number of interacting particles at the beginning of 
% that timestep (interactPrtcls_start_of_t) multiplied by the dynamic lambda2.
%
% It returns an estimate per timestep as well as an overall mean
% calculated from the data beyond 0.2 hours for which the gradient of the
% data is approximately constant.
%
% Units of rates are per hour.

% Calculate tmax_noCC
tmax_noCC = times(end);

% Find when the steady state in interacting particles is reached, if it
% exists
deriv1_interacting = gradient(av_data(1,:),PARAMETERS.tstep_duration);

% Find the timesteps when a constant non-zero gradient in interacting
% particles is reached, if it exists
deriv2_interacting = gradient(deriv1_interacting,PARAMETERS.tstep_duration);
zero_deriv2_interacting_times = times(abs(deriv2_interacting(1:length(times)))<=tol);
zero_deriv2_interacting_times = zero_deriv2_interacting_times...
    (zero_deriv2_interacting_times > 0.2);

% If there are points of constant gradient
if zero_deriv2_interacting_times 
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
    indices = indices -2;
    % Find the average gradient at each timestep. Estimate lambda_2 using 
    %       (gradient of interacting curve) * tstep_duration = 
    %       lambda_1*(# free particles at start of t) - 
    %       lambda_2* (# interacting particles at start of t)
    if length(l1)>1
        l1 = l1(2:end);
    end
    l2_array = (l1 .* freePrtcls_start_of_t(2:end) ...
            - deriv1_interacting(3:end))./...
            interactPrtcls_start_of_t(2:end); % to avoid the initial 0
    %l2_mean = mean(l2_array(indices)); % non-weighted
    [rows,~]=size(av_data);
    if rows>1
        maxT = min(tmax_noCC,(length(l2_array)-1)*PARAMETERS.tstep_duration);
        l2_mean=w8mean(l2_array(1:int64(maxT/PARAMETERS.tstep_duration + 1)),...
            av_data(3,1:int64(maxT/PARAMETERS.tstep_duration + 1))); % weighted
    else
        maxT = min(tmax_noCC,(length(l2_array)-1)*PARAMETERS.tstep_duration);
        l2_mean=w8mean(l2_array(1:int64(maxT/PARAMETERS.tstep_duration + 1)),...
            av_data(1:int64(maxT/PARAMETERS.tstep_duration + 1))); % weighted
    end
else
    fprintf("\nlambda2 was not able to be estimated via differences due to noise \n")
    l2_array = [];
    l2_mean = [];
end
