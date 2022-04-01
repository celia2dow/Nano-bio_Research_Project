function lambdas = input_EWT_from_fraction(fractions,num_hours)
% INPUT_EWT_FROM_FRACTION calculates the L rates that can be input into
% cells_simulation.m given the desired fractions of particles in each of 
% the L stages included in the particle-cell interaction model after so
% many hours. E.g. lambda(1) = rate of association, lambda(2) = rate of 
% internalisation for a model with L=2 stages. The fracitons must be
% numbers between 0 and 1. For L+1 inputs, there are L outputs.
% 
% E.g.      [l1, l2] = input_EWT_from_fraction(...
%                           frac_associated,frac_internalised,num_hours);

lambdas = zeros(1,length(fractions));
guessLAM = 0.1; 
for i = 1:length(fractions)
    frac = fractions(i);
    % Create the n-lambda hypoexponential CDF
    fun_of = hypoexpCDF(i,lambdas,frac);
    fun = @(lambda_i) fun_of(lambda_i,num_hours);
    % Estimate lambda(i) by finding the root of this function (away from 0)
    lambdas(i)=fzero(fun,guessLAM);
end
end

function item = item_j(vec,j)
    % Picks out the jth element of a vector in a single line of code
    item = vec(j);
end
