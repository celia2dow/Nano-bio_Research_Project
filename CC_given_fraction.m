function CC = CC_given_fraction(methods,internalPrtcls_start_of_t,frction, PARAMETERS)
% CC_GIVEN_FRACTION estimates the carrying capacity implemented by solving
% the problem:
%       CC = internalPrtcls_start_of_t/(1-frction)
% where 'frction' depends on which method is being used for estimation:
%       DIFFERENCES: frction = gradient of internalPrtcls/
%                               interactPrtcls_start_of_t
%       DYNAMIC RATE: frction = dynamic lambda 2/ estimated mean lambda 2
%                   # either using calculations via the differences method
%                   # or calculations via the mix method
%
% The CC calculated at each time is returned in addition to a mean estimate
% for the CC calculated from all data beyond 12 hours.

% Estimate CC assuming that the dynamic lambda2 = the original lambda2
% * (1 - #internalised per cell/CC)
rows = length(methods);
for i = 1:rows
    CC_array = internalPrtcls_start_of_t ./ (1-frction(i,:));
    CC_last_12hrs = CC_array(12/PARAMETERS.tstep_duration:end);
    CC_mean = mean(CC_last_12hrs(CC_last_12hrs>0));
    CC.(methods{i}) = CC_array;
    CC.([methods{i} '_mean']) = CC_mean;
end
end