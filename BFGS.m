function [minf, lam_, errCode, itCount, fhist, xhist] = BFGS (m, lam0, preci, maxIt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a_minf, a_lam_, a_errCode, a_itCount, a_fhist, a_xhist] = BFGS([0.5, 0.1], [2, 2], 0.0000001, 1000);
% [a_minf, a_lam_, a_errCode, a_itCount, a_fhist, a_xhist] = BFGS(lam0_true_5, [0, -1, 20, -3, -1], 0.0000001, 3000);
%
% lam0: initial value of lambda 
% preci: precision 
% maxIt: maximum iteration

    %% Initialization
    timer_(-1);
    if size(m, 1) == 1
        m = transpose(m);
    end
    if size(lam0, 1) == 1
        lam0 = transpose(lam0);
    end
    errCode = 0;
    minf = 0;
    lam_ = 0;
    ReadData();
    integral_impl(0, 0, 0, 0);

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
    % phi(alpha) = f(x + alpha * p)
    phi_int = @(alpha, p, x, lam) exp((lam + alpha.*p)'*power(x, 0:n-1)');
    phi = @(alpha, p, lam) integral_impl(@(x) phi_int(alpha, p, x, lam), 0, 1) - (lam+alpha.*p)'*m;
    % derivative of phi(alpha) = f'(x + alpha * p) w.r.t. alpha
    phid_int = @(alpha, p, x, lam) (p'*power(x, 0:n-1)').*exp((lam + alpha.*p)'*power(x, 0:n-1)');
    phid = @(alpha, p, lam) integral_impl(@(x) phid_int(alpha, p, x, lam), 0, 1) - p'*m;

    %% Quasi Newton with BFGS's implementation (Sherman optimized)
    % s_k = x_k - x_{k-1}
    % y_k = F(x_k) - F(x_{k-1})
    %                  sy^T                           ys^T      ss^T     
    % A_k^{-1} = (I - -------) * A_{k-1}^{-1} * (I - ------) + ------
    %                  s^Ty                           s^Ty      s^Ts
    %% First iterations and recorder
    timer_(0);
    lam_last = lam0;
    feval_last = f(lam0);
    y_last = zeros(n, 1);  % last partial (gradient) eval
    for i = 0:n-1
    
        y_last(i+1) = p(i, lam_last);
    end
    % Other initialization
    g_ = zeros(n, n);  % jacobian evaluated at lam0
    for i = 0:n-1
        for j = 0:i
            
            g_(i+1, j+1) = g(i, j, lam_last);
            g_(j+1, i+1) = g_(i+1, j+1);
        end
    end
    %A_ = g_ + (y - y_last - g_*(lam-lam_last)) * (lam-lam_last)' ./ norm(lam-lam_last)^2;
    A_ = g_;
    % Modified Cholosky to ensure SPD
    [L, DMC, P] = modified_cholesky(A_);
    A_ = P'*L*DMC*L'*P;
    A_ = A_^(-1);  % First inverse
    % Manually perform one iteration (using Newtons)
    w = A_ * y_last;
    lam = lam_last - w ./ norm(w);
    feval = f(lam);
    y = zeros(n, 1);  % partial (gradient) eval
    for i = 0:n-1

        y(i+1) = p(i, lam);
    end
    % Recorder and initialization
        % s = zeros(n, 1);  % s_k = x_k - x_{k-1}
        % y_delta = zeros(n, 1);  % delta y
        % w = zeros(n, 1);  % update gap
    it = 1;
    feval_his = zeros(maxIt+1, 1);
    feval_his(1) = feval_last;
    lam_his = zeros(maxIt+1, n);  % x history
    lam_his(1, :) = lam_last;

    %% 2nd to maxIt iterations
    
    % iteration starts
    while it <= maxIt %&& abs(feval - feval_last) > preci

        fprintf('------ iter: %d ------\n', it);
        % Calculate delta x and delta y
        timer_(0);
        s = lam - lam_last;
        y_delta = y - y_last;
        % fprintf('**** Delta y: %f ****\n',norm(y_delta));
        
        % Using BFGS method to update Jacobian
        % A_ is inverse of Jacobian
        A_ = (eye(n)-s*y_delta'/(s'*y_delta)) * A_ * (eye(n)-y_delta*s'/(s'*y_delta))+s*s'/(s'*y_delta);
        if nnz(eig(A_) < 0) 
            fprintf('Not positive definite.\n');
            
            A_ = (A_ + A_') ./ 2;
            [L, DMC, P] = modified_cholesky(A_);
            A_ = P'*L*DMC*L'*P;
        end
        timer_(2);

        w = -A_ * y;
        fprintf('**** Current W Norm: %f ****\n',norm(w));
        [alpha, errCode] = LineSearch(@(alpha) phi(alpha, w, lam), @(alpha) phid(alpha, w, lam), w, it, true);
        if ~(errCode == 0 || errCode == -2), break;end
        w = alpha .* w;
        timer_(3);
        
        % New lam
        lam_his(it+1, :) = lam';
        lam_last = lam;
        lam = lam + w;
        fprintf('**** Delta Lambda Norm: %f ****\n',norm(w));
        timer_(0);

        % re-evaluate f
        feval_last = feval;
        feval_his(it+1) = feval;
        feval = f(lam);
        %if feval < 0 && feval_last > 0
        %    integral_impl(0, 0, 0, 2);
        %end
        if it > 10 && feval_last < 0 && feval > 0  % END

            integral_impl(0, 0, 0, 5);
            it_ = it;
            minf = f(lam_his(it_, :)');
            while minf > 0
                it_ = it_ - 1;
                minf = f(lam_his(it_, :)');
            end
            lam_ = lam_his(it_, :)';
            break;
        end

        timer_(4);

        % re-evaluate y (gradient of f)
        y_last = y;
        for i = 0:n-1

            y(i+1) = p(i, lam);
        end
        timer_(5);

        it = it + 1;
    end

    %% Result
    if it > maxIt, errCode = 1; end
    if isnan(minf), errCode = 2;end
    itCount = it - 1;
    fhist = feval_his(1:it);
    xhist = lam_his(1:it, :);
    timer_(-3);
end
