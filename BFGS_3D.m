function [minf, lam_, errCode, itCount, fhist, xhist] = BFGS_3D (lam0, preci, maxIt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a_minf, a_lam_, a_errCode, a_itCount, a_fhist, a_xhist] = BFGS([0.5, 0.1], [2, 2], 0.0000001, 1000);
% [a_minf, a_lam_, a_errCode, a_itCount, a_fhist, a_xhist] = BFGS(lam0_true_5, [0, -1, 20, -3, -1], 0.0000001, 3000);
%
% lam0: initial value of lambda 
% preci: precision 
% maxIt: maximum iteration

    %% Initialization
    timer_(-1);
    %{
    if size(m, 1) == 1
        m = transpose(m);
    end
    %}
    if size(lam0, 1) == 1
        lam0 = transpose(lam0);
    end
    errCode = 0;
    integral_impl(1, -2);

    n = length(lam0);

    global i1_i2_i3_m;
    data = importdata('data/m3.txt');
    i1_i2_i3_m = data.data;

    %% Quasi Newton with BFGS's implementation (Sherman optimized)
    % s_k = x_k - x_{k-1}
    % y_k = F(x_k) - F(x_{k-1})
    %                  sy^T                           ys^T      ss^T     
    % A_k^{-1} = (I - -------) * A_{k-1}^{-1} * (I - ------) + ------
    %                  s^Ty                           s^Ty      s^Ts
    %% First iterations and recorder
    timer_(0);
    lam_last = lam0;
    feval_last = Gauss_integral3D(lam_last, 'f');
    y_last = Gauss_integral3D(lam_last, 'g');  % last partial (gradient) eval
    
    % Other initialization
    A_ = Gauss_integral3D(lam_last, 'h');
    
    % Modified Cholosky to ensure SPD
    [L, DMC, P] = modified_cholesky(A_);
    A_ = P'*L*DMC*L'*P;
    A_ = A_^(-1);  % First inverse
    % Manually perform one iteration (using Newtons)
    w = A_ * y_last;
    lam = lam_last - w ./ norm(w);
    
    feval = Gauss_integral3D(lam, 'f');
    y = Gauss_integral3D(lam, 'g');  % partial (gradient) eval

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
 
        if it == 492
            fprintf('1');
        end
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
        [alpha, errCode] = LineSearch_3D(lam, w, it, true);
        if errCode ~= 0, break;end
        w = alpha .* w;
        timer_(3);
        
        % New lam
        lam_his(it+1, :) = lam';
        lam_last = lam;
        lam = lam + w;
        fprintf('**** Delta Lambda Norm: %f ****\n',norm(w));
        timer_(0);

        % re-evaluate y (gradient of f)
        y_last = y;
        y = Gauss_integral3D(lam, 'g');
        timer_(4);

        % re-evaluate f
        feval_last = feval;
        feval_his(it+1) = feval;
        feval = Gauss_integral3D(lam, 'f');
        %if feval < -280
        %    integral_impl(1, -1);
        %end

        timer_(5);
        
        it = it + 1;
        
    end
    
    %% Result
    if it > maxIt, errCode = 1; end
    minf = f(lam);
    if isnan(minf), errCode = 2;end
    lam_ = lam;
    itCount = it-1;
    fhist = feval_his(1:it);
    xhist = lam_his(1:it, :);
    timer_(-3);
end
