function [minf, l_1] = Broydens(m, lam0, preci, maxIt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lam0: initial value of lambda 
% preci: precision 
% maxIt: maximum iteration

    %% Initialization
    if size(m, 1) == 1
        m = transpose(m);
    end
    if size(lam0, 1) == 1
        lam0 = transpose(lam0);
    end
    
    n = length(m);
    
    % the integrand part in 'f'
    f_int = @(x, lam, n) exp(lam'*power(x, 0:n-1)');
    % build the core function
    f = @(lam, m, n) integral_impl(@(x) f_int(x, lam, n), 0, 1) - lam'*m;
    % the jacobian function (hessian)
    g_int = @(i, j, x, lam, n) x.^(i+j).*exp(lam'*power(x, 0:n-1)');
    g = @(i, j, lam, n) integral_impl(@(x) g_int(i, j, x, lam, n), 0, 1);
    % partial derivative function (gradient)
    p_int = @(i, x, lam, n) x.^i.*exp(lam'*power(x, 0:n-1)');
    p = @(i, lam, m, n) integral_impl(@(x) p_int(i, x, lam, n), 0, 1) - m(i+1);


    %% Quasi Newton with Broyden's implementation (Sherman optimized)
    % s_k = x_k - x_{k-1}
    % y_k = F(x_k) - F(x_{k-1})
    %                           (s - A^{-1}y) s^T A^{-1}
    % A_k^{-1} = A_{k-1}^{-1} + ------------------------
    %                                s^T A^{-1} y
    % w_k = A_k^{-1} F(x_k)
    % x_{k+1} = x_k - w_k
    
    feval = f(lam0);
    feval_last = abs(feval) / 2;
    lam = ;
    lam_last = lam0;
    y = zeros(n, 1);  % partial (gradient) eval
    y_last = zeros(n, 1);  % last partial (gradient) eval
    s = zeros(n, 1);  % s_k = x_k - x_{k-1}
    % manually perform one iteration
    A = ;
    A_ = ;  % First inverse
    
    it = 1;
    
    feval_his = zeros(maxIt, 1);
    lam_his = zeros(maxIt, n);  % x history
    while it <= maxIt && abs(feval - feval_last) > preci
       
        % Calculate delta x and delta y
        s = lam - lam_last;
        y_delta = y - y_last;
        % Using Sherman Morrison, A_ means A inverse
        A_ = A_ + (s-A_*y_delta) * s'*A_ ./ (s'*A_*y_delta);
        w = A_ * y;
        lam = lam - w;
    end
    
        
    end
    
    
    
    
    
    
end
