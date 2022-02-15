function [lambda1,lambda2] = input_EWT_from_fraction(frac_associated,...
    frac_internalised,num_hours)
% INPUT_EWT_FROM_FRACTION calculates the rates that can be input into
% cells_simulation.m given the desired fractional association
% (frac_associated) and the desired fractional internalisation
% (frac_internalised) after so many hours (num_hours). frac_associated
% and frac_internalised are numbers between 0 and 1.

% Find lambda1 from the CDF of the exponential distribution, F(t), given
% that F(num_hours)=frac_associated
lambda1 = log(1-frac_associated)/(-num_hours); 

% Don't want the root finding function to arrive at lambda1=lambda2
if 1 >= frac_internalised/(frac_associated^2) % lambda1 >= lambda2
    guess = lambda1/3;
else % lambda1 < lambda2
    guess = 3*lambda1;
end

% Create a function from the CDF of the hypoexponential distribution, G(t),
% given that G(num_hours)=frac_internalised and lambda1 is as calculated.
fun = @(lambda2) 1 - frac_internalised - 1/(lambda2-lambda1) * ...
    (lambda2 * exp(-lambda1 * num_hours) - ...
    lambda1 * exp(-lambda2 * num_hours));

% Estimate lambda2 by finding the root of this function (away from lambda1)
lambda2 = fzero(fun,guess);
end
