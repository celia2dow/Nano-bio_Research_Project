function [tmax_noCC, l2_mean, l2_array] = tmax_l2_from_hypoexpCDF(times,...
    hypoCDF,av_data,PARAMETERS,guess,l1_mean,tol_l2)
% TMAX_L2_FROM_HYPOEXPOCDF estimates lambda2 by equating the data (av_data)
% on average numbers of particles internalised over the initial dosage of 
% particles per cell (PARAMETERS.prtcls_per_cell) to the hypoexponential CDF
% (hypoCDF) evaluated at that time (times), given an initial guess (guess)
% that is checked to not give the same rate as the MLE estimate for lambda1
% (l1_mean). The estimates of lambda2 are effective rates that capture how
% the CC reduces the rate of internalisation on each timestep but do not
% equate to the dynamic rate governing successful events on each timestep.
%
% It then calls upon TMAX_FROM_HYPOEXPCDF to estimate the time (tmax_noCC) 
% at which it is assumed that the carrying capacity kicks in significantly
% (as defined by the input tol_l2) if indeed there is a carrying capacity.
%
% It returns tmax_noCC, an estimate of lambda2 per timestep, and an overall 
% mean lambda2 calculated from the data from 1 hour to tmax_noCC. If there 
% is no carrying capacity included, tmax_noCC = PARAMETERS.simulation_duration.
%
% Units of rates are per hour.

% CALCULATE  L2_ARRAY
l2_array = zeros(1,length(times));
smth_intern = smooth(av_data(2,:))';
for i = 1:length(l2_array)
    % Hypoexponential CDF - fraction internalisation
    func = @(l2) hypoCDF(l2,times(i)) - smth_intern(i)/PARAMETERS.prtcls_per_cell;
    % Estimate lambda2 by finding the root of this function (away from lambda1)
    % i.e., by substituting the data into the hypoexponential distribution
    l2_array(i) = fzero(func,guess(2));
    if l2_array(i) == l1_mean
        l2_array(i) = fzero(func,guess(1));
    end
end

% CALCULATE TMAX_NOCC
tmax_noCC = tmax_from_hypoexpCDF(PARAMETERS,l2_array,tol_l2,times);

% CALCULATE L2_MEAN
%l2_mean=mean(l2_array(int64(1/PARAMETERS.tstep_duration + 1):...
%    int64(tmax_noCC/PARAMETERS.tstep_duration + 1))); % non weighted
l2_mean=w8mean(l2_array(1:int64(tmax_noCC/PARAMETERS.tstep_duration + 1)),...
    av_data(3,1:int64(tmax_noCC/PARAMETERS.tstep_duration + 1))); % weighted
end
