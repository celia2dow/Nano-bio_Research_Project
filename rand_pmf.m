function X = rand_pmf(x,PMF,NUM)
% RAND_PMF Random numbers from a user defined discrete distribution using
% the principles of the Inverse CDF transform method of random number 
% generation
%   Adapted from Furlan (2022). Generation of random number given the pmf 
%   (https://www.mathworks.com/matlabcentral/fileexchange/26608-generation-
%   of-random-number-given-the-pmf), MATLAB Central File Exchange. 
%   Retrieved February 15, 2022. 
%
% X=rand_pmf(x,PMF,NUM)
% Input:
% x   : set of the all possible values that the desired random signal can
%       assume
% PMF : vector that cointains the probability of each possible 
%       value of x
% NUM : number of random values to be chosen
% output:
% X   : random signal with the desired pmf
%
% Example: 
% PMF=[1/3 1/3 1/3]
% x=[1 2 3];
% NUM=100;
% X=rand_pmf(x,PMF,NUM);
a=cumsum((PMF(:)))*(1./rand(1,NUM));
b=a>=ones(length(PMF),NUM);
[~,index]=max(b);
x=x(:);
X=x(index);
end
