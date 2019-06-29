function [minf, lam_, errCode, itCount, fhist, xhist] = Broydens(m, lam0, preci, maxIt)
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
    f_int = @(x, lam) exp(lam'*power(x, 0:n-1)');
    % build the core function
    f = @(lam) integral_impl(@(x) f_int(x, lam), 0, 1) - lam'*m;
    % the jacobian function (hessian)
    g_int = @(i, j, x, lam) x.^(i+j).*exp(lam'*power(x, 0:n-1)');
    g = @(i, j, lam) integral_impl(@(x) g_int(i, j, x, lam), 0, 1);
    % partial derivative function (gradient)
    p_int = @(i, x, lam) x.^i.*exp(lam'*power(x, 0:n-1)');
    p = @(i, lam) integral_impl(@(x) p_int(i, x, lam), 0, 1) - m(i+1);


    %% Quasi Newton with Broyden's implementation (Sherman optimized)
    % s_k = x_k - x_{k-1}
    % y_k = F(x_k) - F(x_{k-1})
    %                           (s - A^{-1}y) s^T A^{-1}
    % A_k^{-1} = A_{k-1}^{-1} + ------------------------
    %                                s^T A^{-1} y
    % w_k = A_k^{-1} F(x_k)
    % x_{k+1} = x_k - w_k
    
    lam_last = lam0;
    feval_last = f(lam0);
    y_last = zeros(n, 1);  % last partial (gradient) eval
    for i = 0:n-1
    
        y_last(i+1) = p(i, lam_last);
    end
    % manually perform one iteration
    lam = lam_last + 2.*rand(n, 1);
    feval = f(lam);
    y = zeros(n, 1);  % partial (gradient) eval
    for i = 0:n-1

        y(i+1) = p(i, lam);
    end
    % other initialization
    g_ = zeros(n, n);  % jacobian evaluated at lam0
    for i = 0:n-1
        for j = 0:n-1
            
            g_(i+1, j+1) = g(i, j, lam_last);
        end
    end
    A_ = g_ + (y - y_last - g_*(lam-lam_last)) * (lam-lam_last)' ./ norm(lam-lam_last)^2;
    A_ = A_^(-1);  % First inverse
    s = zeros(n, 1);  % s_k = x_k - x_{k-1}
    y_delta = zeros(n, 1);  % delta y
    w = zeros(n, 1);  % update gap
    it = 1;
    feval_his = zeros(maxIt, 1);
    lam_his = zeros(maxIt, n);  % x history
    
    % iteration starts
    while it <= maxIt && abs(feval - feval_last) > preci
       
        % Calculate delta x and delta y
        s = lam - lam_last;
        y_delta = y - y_last;
        % Using Sherman Morrison, A_ means A inverse
        A_ = A_ + (s-A_*y_delta) * s'*A_ ./ (s'*A_*y_delta);
        w = A_ * y;
        
        % New lam
        lam_his(it, :) = lam';
        lam = lam - w;
        
        % re-evaluate
        feval_last = feval;
        feval_his(it) = feval;
        feval = f(lam);
        
        it = it + 1;
    end
    
    %% Result
    errCode = 0;
    if it > maxIt
        errCode = 1; 
    end
    minf = f(lam);
    if isnan(minf)
        errCode = 2;
    end
    lam_ = lam;
    itCount = it-1;
    fhist = feval_his(1:it-1);
    xhist = lam_his(1:it-1, :);
end
