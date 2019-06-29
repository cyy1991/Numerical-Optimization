function [minf, lam_, errCode, itCount, fhist, xhist] = Newtons (m, lam0, preci, maxIt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Example
% m:       [0.5, 1, 2, 0.7] both row or column vector is fine
% lam0:    [1, 1, 1, 2]
% preci:   terminate condition '|diff of evaluations| < preci'
% maxIt:   max number of iterations
%
% minf:    the resultant minimum function value
% lam_:    the lam values to achieve the minimum value
% errCode: 0 normal exit; 1 exceeds maxIt; 2 fatal error
% itCount: iterations used
% fhist:   history values of f [it * 1]
% xhist:   history choices of x (lambda) [it * n]

    %% Initialization
    timer_(-1);
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

    %% Newton's Method
    %  y = F(x)
    %  w = J^{-1} y
    %  x = x + w
    it = 1;
    feval = f(lam0);
    feval_last = abs(feval) / 2;
    if abs(feval - feval_last) <= preci
        feval_last = 100;  % in case first eval is small
    end
    lam = lam0;
    J = zeros(n, n);
    y = zeros(n, 1);
    feval_his = zeros(maxIt, 1);
    lam_his = zeros(maxIt, n);
    while it <= maxIt && abs(feval - feval_last) > preci

        % Calculate Jacobian
        timer_(0);
        for i = 0:n-1
            for j = 0:n-1
            
                J(i+1, j+1) = g(i, j, lam);
                J(j+1, i+1) = J(i+1, j+1);
            end
            timer_(1);
            y(i+1) = -p(i, lam);  % minus symbol to enable l = l + w
            timer_(2);
        end
        
        % Modified Cholosky to ensure SPD
        [L, DMC, P] = modified_cholesky(J);
        J = P'*L*DMC*L'*P;
        timer_(3);
        
        % Calculate Inverse/Solve Guassian
          % J_inv = inverse_impl(J);
        w = linsolve_impl(J, y);
        timer_(4);
          
        % New lam
        lam_his(it, :) = lam';
          % lam = lam - J_inv * lam;
        lam = lam + w;

        % re-evaluate
        feval_last = feval;
        feval_his(it) = feval;
        feval = f(lam);
        timer_(5);

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
    timer_(-3);
end
