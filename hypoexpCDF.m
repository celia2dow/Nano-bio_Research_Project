function fun = hypoexpCDF(n,lambdas,frac)
% HYPOEXPCDF creates a function hangle for the n-lambda hypoexponential CDF 
% with lambda_n being the only unknown rate. The other known rates are
% contained in lambdas and frac is the fraction between 0 and 1 that the
% CDF takes (i.e. F(t)).
    if n == 1
        % Find lambda(1) from the CDF of the exponential distribution, 
        % F(t), given that F(t)=frac
        fun = @(lambda_n,t) 1 - frac - exp(-lambda_n*t);
    else
        % Create a function from the CDF of the hypoexponential 
        % distribution with n rates, G(t), given that G(t)=frac
        % and lambda(1), ..., lambda(n-1) are contained in lambdas
        fun = @(lambda_n,t) 1 - frac;
        for j = 1:n
            % Weighting for -exp(-lambda(j)*t) is the product of
            % lambda(k) for k ~= j, divided by the product of 
            % (lambda(k)-lambda(j)) for k ~= j.
            weighting = @(lambda_n) prod(lambdas(1:n-1)) * ...
            lambda_n/item_j([lambdas(1:n-1),lambda_n],j);
            for k = [1:j-1 j+1:n]
                weighting = @(lambda_n) weighting(lambda_n) / ...
                    (item_j([lambdas(1:n-1),lambda_n],k)-...
                    item_j([lambdas(1:n-1),lambda_n],j));
            end
            fun = @(lambda_n,t) fun(lambda_n,t) - weighting(lambda_n) * ...
                exp(-item_j([lambdas(1:n-1),lambda_n],j)*t);
        end
    end
end 
