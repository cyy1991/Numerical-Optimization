function [minf, l_1] = BFGS(m, lam0, preci, maxIt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f: target function
% grad: gradient of f
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
    % initial two steps created by simple gradient descent
    n=length(lam0);
    y_2 = zeros(n, 1); % y_k-2
    y_1 = zeros(n, 1); % y_k-1
    for i = 0:n-1

        y_2(i+1) = grad(i, lam0, m, n) ;
    end
    y_2(y_2 > 100) = 100;
    y_2(y_2 < -100) = -100;
    l_2 = lam0 - y_2 * 0.01; % l_k-2
    for i = 0:n-1
        
        y_2(i+1) = grad(i, l_2, m, n) ;
    end
    y_2(y_2 > 100) = 100;
    y_2(y_2 < -100) = -100;
    l_1 = l_2 - y_2 * 0.01; % l_k-1
    for i = 0:n-1
        
        y_1(i+1) = grad(i, l_1, m, n) ;
    end
    y_1(y_1 > 90) = 90;
    y_1(y_1 < -90) = -90;

    f_last = f(l_2, m, n);
    f_cur = f(l_1, m, n);
    % last step hessian  (could be replaced?)
    D_last = eye(n);  % D denotes the inverse of hessian

    %% Quasi Newton
    iter = 0;
    delta_y = y_1 - y_2;
    while abs(f_cur - f_last) > preci  % note precision is used when comparing gradient change to save computation cost of f(l), may not be optimal 
        
        f_last = f_cur;
        s = l_1 - l_2;
        D_last = D_last + (s - D_last * delta_y) * s' * D_last / (s' * D_last * delta_y);
        %D_last = D_last + (delta_y - D_last * s)/ (s'*s) * s';
        w =  D_last * y_1;
        l_2 = l_1;
        l_1 = l_1 - 0.05 * w;
        y_2 = y_1;
        for i = 0:n-1
            
            y_1(i+1) = grad(i, l_1, m, n) ;
        end
        iter = iter + 1;
        if iter > maxIt
            
            warning('reach maximum iteration, process terminated');
            break 
        end
        delta_y = y_1 - y_2;
        f_cur = f(l_1, m, n);
    end

    display(iter)
    minf = f(l_1, m, n);
end
