function intern_CC = analytic_distrib_with_CC(times,lam,CC,PARAMETERS, ...
    hypoCDF,expCDF,assocPrtcls,L)
% ANALYTIC_DISTRIB_WITH_CC creates an array of the analytical distribution
% of the internalised particles (intern_curve) when given heuristic 
% estimates for lambdas and the arrying capacity (CC). 
intern_CC = zeros(1,length(times));
new_interns_CC = zeros(1,length(times));
if L == 1
    % Assume that new internalisations are the product of the number of
    % free particles available multiplied by the exponential distribution
    % evaluated at the length of a timestep.
    new_interns_CC(1) = PARAMETERS.prtcls_per_cell .* expCDF(lam,0);
    intern_CC(1) = new_interns_CC(1);
    for i = 2:length(times)
        l1_dyn = lam*(1-intern_CC(i-1)/CC);
        freePrtcls = PARAMETERS.prtcls_per_cell - intern_CC(i-1);
        new_interns_CC(i) = expCDF(l1_dyn,times(2)) .* freePrtcls;
        intern_CC(i) = intern_CC(i-1) + new_interns_CC(i);
    end
elseif L == 2
    % It is assumed that the new internalisation within a timestep will be a 
    % combination of new internalisation events from free (described by the 
    % hypoexponential distirbution, hypoCDF, scaled by the remaining dosage) 
    % and from bound (described by the exponential distribution, expCDF, scaled
    % by the number bound).
    new_interns_CC(1) = PARAMETERS.prtcls_per_site .* hypoCDF(lam,0);
    intern_CC(1) = new_interns_CC(1);
    for i = 2:length(times)
        l2_dyn = lam*(1-intern_CC(i-1)/CC);
        freePrtcls = PARAMETERS.prtcls_per_site - assocPrtcls(i-1);
        boundPrtcls = assocPrtcls(i-1)- intern_CC(i-1);
        new_interns_CC(i) = hypoCDF(l2_dyn,times(2)) .* freePrtcls +...
            expCDF(l2_dyn,times(2)) .* boundPrtcls;
        intern_CC(i) = intern_CC(i-1) + new_interns_CC(i);
    end
end
end