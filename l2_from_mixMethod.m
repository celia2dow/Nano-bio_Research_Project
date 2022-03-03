function [l2_mean, l2_array] = l2_from_mixMethod(times,hypoCDF,expCDF,...
    av_data,freePrtcls_start_of_t, interactPrtcls_start_of_t,...
    PARAMETERS,guess,l1_mean,tmax_noCC)
% L2_FROM_MIXMETHOD estimates lambda2 by equating the data (av_data) on 
% average number of new particles internalised at a time (times) to the
% free particles available at the beginning of the timestep 
% (freePrtcls_start_of_t) multiplied by the hypoexponential CDF (hypoCDF) 
% evaluated at PARAMETERS.tstep_duration plus the number of particles 
% interacting at the beginning of the timestep (interactPrtcls_start_of_t) 
% multiplied by the exponential CDF (expCDF) of rate lambda2 evaluated at 
% PARAMETERS.tstep_duration. An initial guess (guess) is given that is
% checked to not give the same rate as the MLE estimate for lambda1
% (l1_mean). The estimates of lambda2 are dynamic rates that equate to the
% reduced rate governing successful events on each timestep. 
%
% It returns an estimate per timestep as well as an overall mean
% calculated from the data from 1 hour to tmax_noCC hours (given the
% duration of a timestep PARAMETERS.tstep_duration).
%
% Units of rates are per hour.

% CALCULATE L2_ARRAY
l2_array = zeros(1,length(times)-1);
new_interns = av_data(2,2:end) - av_data(2,1:(end-1));
for i = 1:length(l2_array)
    % F(tstep_duration) * freePrtcls_start_of_t + G(tstep_duration) *
    % interactPrtcls_start_of_t - new_interns
    %   F(t) = Hypoexponential CDF with rates l1 and dynamic l2
    %   G(t) = Exponential CDF rates with rate dynamic l2
    func = @(l2) hypoCDF(l2)  * freePrtcls_start_of_t(i) + ...
        expCDF(l2) * interactPrtcls_start_of_t(i) - new_interns(i);
    % Estimate lambda2 by finding the root of this function (away from lambda1)
    % i.e., by substituting the data into the hypoexponential distribution
    l2_array(i) = fzero(func,guess(2));
    if l2_array(i) == l1_mean
        l2_array(i) = fzero(func,guess(1));
    end
end

% CALCULATE L2_MEAN
%l2_mean=mean(l2_array(int64(1/PARAMETERS.tstep_duration + 1):...
%    int64(tmax_noCC/PARAMETERS.tstep_duration + 1))); % non weighted
maxT = min(tmax_noCC,(length(l2_array)-1)*PARAMETERS.tstep_duration);
l2_mean=w8mean(l2_array(1:int64(maxT/PARAMETERS.tstep_duration + 1)),...
    av_data(3,1:int64(maxT/PARAMETERS.tstep_duration + 1))); % weighted
end
