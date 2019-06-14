function [minf, lam_] = Newtons (m, lam0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Example
% m:     [0.5, 1, 2, 0.7]
% lam0:  [1, 1, 1, 1]
%
% minf:  the resultant minimum function value
% lam_:  the lam values to achieve the minimum value
%

    n = length(m);

    % the integrand part in 'f'
    f_int = @(x, lam) exp(lam*power(x, 0:n-1)');
    % build the core function
    f = @(lam) integral_impl(@(x) f_int(x, lam), 0, 1) - lam*m';
    
    minf = f(lam0);
    lam_ = lam0;
end
