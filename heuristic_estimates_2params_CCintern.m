% SCRIPT for making heuristic estimates in the 2 parameter case with or
% without a carrying capacity on the last transition (internalisation).

% HEURISTIC ESTIMATES
%   of tmax_noCC, lambda1, lambda2 and CC
% Use the provided data on the mean numbers of particles per cell that are 
% free, interacting and internalised and estimate the parameters of the 
% exponential and hypoexponential distributions that the inter-event times 
% follow (of particle association and particle internalisation respectively).

fprintf("\nHeuristic estimates for lambda from free to stage 1 and from stage 1 to stage 2: \n(units are per hour)\n")

% The mean number of free particles available to a cell at the beginning of
% each timestep is the number that haven't been associated on the previous
% timestep divided amongst the lattice sites. 
freePrtcls_start_of_t = (PARAMETERS.prtcls_per_cell * PARAMETERS.culture_dim^2 .* ...
    ones(1,length(binrng)) - (means(3,:) .* total.cell_population) ./ num_runs) ./ ...
    (PARAMETERS.culture_dim^2);
% Free particles left over when association is the mean data minus 1
% standard deviation (std. dev.) - i.e. lower bound association and thus
% upper bound free particles
freePrtcls_start_of_tUP = (PARAMETERS.prtcls_per_cell * PARAMETERS.culture_dim^2 .* ...
    ones(1,length(binrng)) - (lower1(3,:) .* total.cell_population) ./ num_runs) ./ ...
    (PARAMETERS.culture_dim^2);
% Free particles left over when association is the mean data plus 1
% standard deviation (std. dev.) - i.e. upper bound association and thus
% lower bound free particles.
freePrtcls_start_of_tLO = (PARAMETERS.prtcls_per_cell * PARAMETERS.culture_dim^2 .* ...
    ones(1,length(binrng)) - (upper1(3,:) .* total.cell_population) ./ num_runs) ./ ...
    (PARAMETERS.culture_dim^2);
% No interaction occur prior to 0 hours, so the entry for '0 hours' 
% (timestep 0 or typically index 1) is in fact the free particles available 
% at the beginning of timestep 1 (typically index 2) from 0-0.1667 hours 
% and so on. The last entry is simply the number of free particles 
% available per cell at the end of the simulation.
freePrtcls_start_of_t = freePrtcls_start_of_t(1:(end-1)); % for timestep 1, 2, ...
freePrtcls_start_of_tUP = freePrtcls_start_of_tUP(1:(end-1)); 
freePrtcls_start_of_tLO = freePrtcls_start_of_tLO(1:(end-1)); 

% The mean number of particles interacting with a cell at the beginning of
% each timestep is the number that are recorded in the previous timestep.
% I.e. the interacting particles available at the beginning of timestep 1
% (typically index 2) from 0-0.1667 hours is the number recorded as
% interacting at the end of timestep 0 (typically index 1) at 0 hours.
interactPrtcls_start_of_t = means(1,1:(end-1)); % for timestep 1, 2, ...
% Interacting particles corresponding to upper bound association and lower
% bound internalisation - i.e. upper bound interacting particles
interactPrtcls_start_of_tUP = upper1(3,:)-lower1(2,:); % full array
% Interacting particles corresponding to lower bound association and upper
% bound internalisation - i.e. lower bound interacting particles
interactPrtcls_start_of_tLO = lower1(3,:)-upper1(2,:); % full array

% ESTIMATE LAMBDA 1
%   via DIFFERENCES method
% Estimate lambda_1, the parameter for the exponential distribution that
% describes the inter-association times of particles per cell, by
% calculating the mean fraction of association events that occur (the
% number that do occur in a timestep over the number that could occur -
% i.e., the number of free particles) over the duration of time passed.
est_lambda1.using_diffs = (means(3,2:end)- means(3,1:(end-1)))./... new in tstep 1, 2, ...
    (binrng(2).*freePrtcls_start_of_t); % free at the beginning of tstep 1, 2, ...
est_lambda1.using_diffs_mean = w8mean(est_lambda1.using_diffs,means(3,2:end));
% Upper bound
est_lambda1UP.using_diffs = (upper1(3,2:end)- upper1(3,1:(end-1)))./... 
    (binrng(2).*freePrtcls_start_of_tLO); 
est_lambda1UP.using_diffs_mean = w8mean(est_lambda1UP.using_diffs,upper1(3,2:end));
% Lower bound
est_lambda1LO.using_diffs = (lower1(3,2:end)- lower1(3,1:(end-1)))./... 
    (binrng(2).*freePrtcls_start_of_tUP); 
est_lambda1LO.using_diffs_mean = w8mean(est_lambda1LO.using_diffs,lower1(3,2:end));

% ESTIMATE LAMBDA 1
%   via MLE method
smth_assoc = smooth(means(3,:))'; % number associated at end of tstep 0, 1, 2, ...
est_lambda1.MLE = smth_assoc(2:end)./ ... % num assoc at end of tstep 1, 2, ...
    (binrng(2:end).*PARAMETERS.prtcls_per_cell);% time at end of tstep 1, 2, ...
est_lambda1.MLE_mean = w8mean(est_lambda1.MLE,means(3,2:end));
% Upper bound
smth_assocUP = smooth(upper1(3,:))'; 
est_lambda1UP.MLE = smth_assocUP(2:end)./ ... 
    (binrng(2:end).*PARAMETERS.prtcls_per_cell);
est_lambda1UP.MLE_mean = w8mean(est_lambda1UP.MLE,upper1(3,2:end));
% Lower bound
smth_assocLO = smooth(lower1(3,:))'; 
est_lambda1LO.MLE = smth_assocLO(2:end)./ ...
    (binrng(2:end).*PARAMETERS.prtcls_per_cell); 
est_lambda1LO.MLE_mean = w8mean(est_lambda1LO.MLE,lower1(3,2:end));

% CDF of hypoexponential distribution
CDF.hypoexp = @(l1, l2, t) 1 - 1./(l2-l1) .* (l2 .* exp(-l1 .* t) - l1 .* exp(-l2 .* t));
CDF.hypoexp_MLEl1 = @(l2,t) CDF.hypoexp(est_lambda1.MLE_mean,l2,t); 
CDF.hypoexp_MLEl1UP = @(l2,t) CDF.hypoexp(est_lambda1UP.MLE_mean,l2,t); % Upper bound
CDF.hypoexp_MLEl1LO = @(l2,t) CDF.hypoexp(est_lambda1LO.MLE_mean,l2,t); % Lower bound
CDF.hypoexp_MLEl1_tstep = @(l2) CDF.hypoexp_MLEl1(l2,PARAMETERS.tstep_duration); 
CDF.hypoexp_MLEl1UP_tstep = @(l2) CDF.hypoexp_MLEl1UP(l2,PARAMETERS.tstep_duration); % Upper bound
CDF.hypoexp_MLEl1LO_tstep = @(l2) CDF.hypoexp_MLEl1LO(l2,PARAMETERS.tstep_duration); % Lower bound
% CDF of exponential distribution
CDF.exp = @(l, t) 1 - exp(- l .*t);
CDF.exp_tstep = @(l) CDF.exp(l,PARAMETERS.tstep_duration);
% Two guesses to feed into fzero in case one gives lambda2=lambda1
guess = [est_lambda1.MLE_mean/3 est_lambda1.MLE_mean*3]; % mean
guessUP = [est_lambda1UP.MLE_mean/3 est_lambda1UP.MLE_mean*3]; % Upper bound
guessLO = [est_lambda1LO.MLE_mean/3 est_lambda1LO.MLE_mean*3]; % Lower bound

% ESTIMATE TIME AT WHICH CC KICKS IN (TMAX_NOCC) and LAMBDA 2
%   via DISTRIBUTION method
[tmax_noCC,est_lambda2.distrib_mean,est_lambda2.distrib] = ... 
    tmax_l2_from_hypoexpCDF(binrng,CDF.hypoexp_MLEl1,means,...
    PARAMETERS,guess,est_lambda1.MLE_mean,tol_l2);
% Upper bound lambda2 (assuming lower bound lambda1)
[~,est_lambda2UP.distrib_mean,est_lambda2UP.distrib] = ... 
    tmax_l2_from_hypoexpCDF(binrng,CDF.hypoexp_MLEl1LO,upper1,...
    PARAMETERS,guessLO,est_lambda1LO.MLE_mean,tol_l2);
% Lower bound lambda2 (assuming upper bound lambda1)
[~,est_lambda2LO.distrib_mean,est_lambda2LO.distrib] = ... 
    tmax_l2_from_hypoexpCDF(binrng,CDF.hypoexp_MLEl1UP,lower1,...
    PARAMETERS,guessUP,est_lambda1UP.MLE_mean,tol_l2);

% ESTIMATE LAMBDA 2
%   via MIX method
[est_lambda2.mix_mean,est_lambda2.mix] = l2_from_mixMethod(...
    binrng,CDF.hypoexp_MLEl1_tstep,CDF.exp_tstep,...
    means,freePrtcls_start_of_t,interactPrtcls_start_of_t,...
    PARAMETERS,guess,est_lambda1.MLE_mean,tmax_noCC);
% Upper bound lambda2 (assuming lower bound lambda1)
[est_lambda2UP.mix_mean,est_lambda2UP.mix] = l2_from_mixMethod(...
    binrng,CDF.hypoexp_MLEl1LO_tstep,CDF.exp_tstep,...
    upper1,freePrtcls_start_of_tUP,interactPrtcls_start_of_tLO(1:end-1),...
    PARAMETERS,guessLO,est_lambda1LO.MLE_mean,tmax_noCC);
% Lower bound lambda2 (assuming upper bound lambda1)
[est_lambda2LO.mix_mean,est_lambda2LO.mix] = l2_from_mixMethod(...
    binrng,CDF.hypoexp_MLEl1UP_tstep,CDF.exp_tstep,...
    lower1,freePrtcls_start_of_tLO,interactPrtcls_start_of_tUP(1:end-1),...
    PARAMETERS,guessUP,est_lambda1UP.MLE_mean,tmax_noCC);

tsteps = 0:PARAMETERS.tstep_duration:tmax_noCC; 

% ESTIMATE LAMBDA 2
%   via DIFFERENCES method
[est_lambda2.using_diffs_mean,est_lambda2.using_diffs] = ...
    l2_from_diffsMethod(tsteps,means,freePrtcls_start_of_t, ...
    interactPrtcls_start_of_t,PARAMETERS,est_lambda1.using_diffs_mean);
% Upper bound lambda2 (assuming lower bound lambda1)
[est_lambda2UP.using_diffs_mean,est_lambda2UP.using_diffs] = ...
    l2_from_diffsMethod(tsteps,interactPrtcls_start_of_tLO,freePrtcls_start_of_tUP, ...
    interactPrtcls_start_of_tLO(1:end-1),PARAMETERS,est_lambda1LO.using_diffs);
% Lower bound lambda2 (assuming upper bound lambda1)
[est_lambda2LO.using_diffs_mean,est_lambda2LO.using_diffs] = ...
    l2_from_diffsMethod(tsteps,interactPrtcls_start_of_tUP,freePrtcls_start_of_tLO, ...
    interactPrtcls_start_of_tUP(1:end-1),PARAMETERS,est_lambda1UP.using_diffs);

% CALCULATE ACTUAL RATES
if PARAMETERS.EWTs_internalise.input_type == "EWT"
    % A list of L rates of particles transitioning between stages of the cell-
    % particle interaction model unaffected by carrying capacity (e.g. free to 
    % stage 1, stage 1 to stage 2, ..., stage L-1 to stage L or internalised).
    % Can be 1 rate per transition or can be 1 rate per cell phase per
    % transition - columns inicate interaction stage and rows indicate cell
    % phase.
    l1 = 1./PARAMETERS.EWTs_internalise.values(1);
    l2 = 1./PARAMETERS.EWTs_internalise.values(2);
elseif PARAMETERS.EWTs_internalise.input_type == "fraction"
    % Calculate the rates from the desired observed fraction
    % associated/internalised at confluence without carrying capacity.
    [l1,l2] = input_EWT_from_fraction(...
        PARAMETERS.EWTs_internalise.values(1),...
        PARAMETERS.EWTs_internalise.values(2),...
        PARAMETERS.EWTs_internalise.values(3));
elseif PARAMETERS.EWTs_internalise.input_type == "prob_and_rates"
    % If the input are just the rates (per hour) from one stage to the
    % next, with the first value being the probability of binding once hit
    % (and therefore not needing scaling by the timestep duration)
    l1 = rate_diffus * PARAMETERS.EWTs_internalise.values(1);
    l2 = PARAMETERS.EWTs_internalise.values(2);
end


% PRINT RATES
fprintf("\nLAMBDA 1: \nACTUAL %5.4e",l1)
fprintf('\nMETHOD \t\tMEAN \t\tLOWER ESTIMATE \tUPPER ESTIMATE')
fprintf('\nDifferences \t%5.4e \t%5.4e \t%5.4e',...
    est_lambda1.using_diffs_mean, est_lambda1LO.using_diffs_mean, ...
    est_lambda1UP.using_diffs_mean);
fprintf('\nMLE \t\t%5.4e \t%5.4e \t%5.4e \n',...
    est_lambda1.MLE_mean, est_lambda1LO.MLE_mean, ...
    est_lambda1UP.MLE_mean);

fprintf("\nLAMBDA 2: \nACTUAL %5.4e",l2)
fprintf('\nMETHOD \t\tMEAN \t\tLOWER ESTIMATE \tUPPER ESTIMATE')
fprintf('\nDifferences \t%5.4e \t%5.4e \t%5.4e',...
    est_lambda2.using_diffs_mean, est_lambda2LO.using_diffs_mean, ...
    est_lambda2UP.using_diffs_mean);
fprintf('\nDistribution \t%5.4e \t%5.4e \t%5.4e',...
    est_lambda2.distrib_mean, est_lambda2LO.distrib_mean, ...
    est_lambda2UP.distrib_mean);
fprintf('\nMix \t\t%5.4e \t%5.4e \t%5.4e\n',...
    est_lambda2.mix_mean, est_lambda2LO.mix_mean, ...
    est_lambda2UP.mix_mean);

% ESTIMATE CARRYING CAPACITY (CC) **if there is one**
%   via DIFFERENCES method
%   via DYNAMIC RATE calculated from DIFFERENCES METHOD
%   via DYNAMIC RATE calculated from MIX METHOD
if any(PARAMETERS.max_prtcls ~= inf) 
    % The number of internalised particles per cell at the start of a 
    % timestep is the number of internalised particles per cell at the end
    % of the previous timestep.
    internalPrtcls_start_of_t = means(2,1:(end-1)); % for timestep 1, 2, ...
    internalPrtcls_start_of_tUP = upper1(2,1:(end-1)); % for timestep 1, 2, ...
    internalPrtcls_start_of_tLO = lower1(2,1:(end-1)); % for timestep 1, 2, ...
    
    % Prepare input for CC_GIVEN_FRACTION
    using_diffs = gradient(means(2,2:end),PARAMETERS.tstep_duration) ./ ... 
        (interactPrtcls_start_of_t .* est_lambda2.using_diffs_mean);
    dynamic_diffs= [0 est_lambda2.using_diffs]./est_lambda2.using_diffs_mean;
    dynamic_mix= est_lambda2.mix./est_lambda2.mix_mean;
    frction = [using_diffs; dynamic_diffs; dynamic_mix];
    methods = {'using_diffs','dynamic_diffs','dynamic_mix'};

    est_CC = CC_given_fraction(methods,internalPrtcls_start_of_t,frction, PARAMETERS);

    % For the upper bound
    using_diffsUP = gradient(upper1(2,2:end),PARAMETERS.tstep_duration) ./ ... 
        (interactPrtcls_start_of_tLO(2:end) .* est_lambda2UP.using_diffs_mean);
    dynamic_diffsUP= [0 est_lambda2UP.using_diffs]./est_lambda2UP.using_diffs_mean;
    dynamic_mixUP= est_lambda2UP.mix./est_lambda2UP.mix_mean;
    frctionUP = [using_diffs; dynamic_diffs; dynamic_mix];
    est_CCUP = CC_given_fraction(methods,internalPrtcls_start_of_tUP,frctionUP, PARAMETERS);

    % For the lower bound
    using_diffsLO = gradient(lower1(2,2:end),PARAMETERS.tstep_duration) ./ ... 
        (interactPrtcls_start_of_tUP(2:end) .* est_lambda2LO.using_diffs_mean);
    dynamic_diffsLO= [0 est_lambda2UP.using_diffs]./est_lambda2UP.using_diffs_mean;
    dynamic_mixLO= est_lambda2UP.mix./est_lambda2UP.mix_mean;
    frctionLO = [using_diffs; dynamic_diffs; dynamic_mix];
    est_CCLO = CC_given_fraction(methods,internalPrtcls_start_of_tLO,frctionLO, PARAMETERS);
    
    % PRINT CARRYING CAPACITIES
    fprintf("\n\nHeuristic estimates for carrying capacity (CC): \n(units are particles)\n")
    
    fprintf("\nCARRYING CAPACITY: \nACTUAL %5.3f",PARAMETERS.max_prtcls(end))
    fprintf('\nMETHOD \t\t\t\tMEAN \tLOWER ESTIMATE \tUPPER ESTIMATE')
    fprintf('\nDifferences \t\t\t%5.3f \t%5.3f \t\t%5.3f',...
        est_CC.using_diffs_mean, est_CCLO.using_diffs_mean, est_CCUP.using_diffs_mean);
    fprintf('\nDynamic using differences \t%5.3f \t%5.3f \t\t%5.3f',...
        est_CC.dynamic_diffs_mean, est_CCLO.dynamic_diffs_mean, est_CCUP.dynamic_diffs_mean);
    fprintf('\nDynamic using mix \t\t%5.3f \t%5.3f \t\t%5.3f\n',...
        est_CC.dynamic_mix_mean, est_CCLO.dynamic_mix_mean, est_CCUP.dynamic_mix_mean);
end
