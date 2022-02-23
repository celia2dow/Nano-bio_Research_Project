function tmax_noCC = tmax_from_hypoexpCDF(PARAMETERS,l2_array,tol_l2,times)
% TMAX_FROM_HYPOEXPCDF estimates the time (tmax_noCC) at which it is
% assumed that the carrying capacity kicks in significantly (as defined by
% the input tol_l2) if indeed there is a carrying capacity. It searches the
% effective lambda2 rates as calculated using the hypoexponential
% distribution (l2_array) and finds the time (times) at which the gradient
% drops by more than the tolerance (tol_l2).

% If there is a carrying capacity invoked, and only assuming that it is a
% carrying capacity on the last stage of the cell-particle interaction
% model.
if any(PARAMETERS.max_prtcls ~= inf)     
    % Assume that the CC begins to have a significant impact on the data
    % when the lambda2 value estimated from the hypoexponential 
    % distribution deviates away from constant by a specific tolerance
    deriv1_l2_distrib = gradient(l2_array);
    bool_slight_drop = zeros(1,length(deriv1_l2_distrib));
    for i = 1:length(deriv1_l2_distrib)
        bool_slight_drop(i) = abs(deriv1_l2_distrib(i))<tol_l2 && ... 
            deriv1_l2_distrib(i)<0;
    end
    indices_slight_drop = find(bool_slight_drop);
    tmax_noCC = times(min(indices_slight_drop(indices_slight_drop>(3/PARAMETERS.tstep_duration))));
    tmax_noCC = max([tmax_noCC, 4]); % Just in case it's too small 
% If there is no carrying capacity
else
    fprintf('No carrying capacity is implemented \n')
    tmax_noCC = PARAMETERS.simulation_duration;
end
end