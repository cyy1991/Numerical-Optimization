function [alpha_star, errCode] = LineSearch_3D(lam, w, it, verbose)

    if nargin<3, verbose = false;end
    c1 = 0.0001;
    c2 = 0.9;

    maxIt = 101;
    multiplier = 1.1482;
    
    n = length(lam);
    
    alpha = zeros(maxIt+1, 1);
    lambda = zeros(n, maxIt+1);
    alpha(1) = 0;
    alpha(2) = 1;
    lambda(:,1) = lam + alpha(1).*w;
    lambda(:,2) = lam + alpha(2).*w;
    % alpha_max = multiplier^(maxIt-1)
    errCode = 0;

    phi_0 = Gauss_integral3D(lambda(:,1), 'f');
    phid_0 = Gauss_integral3D(lambda(:,1), 'g')' * w;
    phi_last = 0;

    %if verbose, fprintf('------ iter: %d ------\n', it);end
    temp = Gauss_integral3D(lambda(:,2), 'f');
    while isnan(temp) || (temp > 10e100)
        
        alpha(2) = alpha(2) / 10;
        lambda(:,2) = lam + alpha(2).*w;
        temp = Gauss_integral3D(lambda(:,2), 'f');
        if alpha(2) < 10e-10
            alpha(2) = 1;
            break;
        end
        if verbose, fprintf('{>>}');end
    end
    while isnan(temp) || (temp > 10e100)
        
        alpha(2) = alpha(2) + 0.05;
        lambda(:,2) = lam + alpha(2).*w;
        temp = Gauss_integral3D(lambda(:,2), 'f');
        if verbose, fprintf('{<<}');end
    end
    for i = 2:maxIt
    
        % evaluate phi(alpha)
        lambda(:,i) = lam + alpha(i).*w;
        phi_a = Gauss_integral3D(lambda(:,i), 'f');
        if verbose, fprintf('alpha:%f  phi:%f\n', alpha(i), phi_a);end
        if (phi_a > phi_0 + c1 * alpha(i) * phid_0) || (i > 2 && phi_a >= phi_last)

            [alpha_star, errCode] = zoom(lam, w, alpha(i-1), alpha(i), verbose);
            if verbose
                lambda_star = lam + alpha_star.* w;
                phi_star = Gauss_integral3D(lambda_star, 'f');
                fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi_star); 
            end
            return;
        end
        % evaluate phid(alpha)
        phid_a = Gauss_integral3D(lambda(:,i), 'g')' * w;
        if abs(phid_a) <= -c2 * phid_0
        
            alpha_star = alpha(i);
            if verbose
                lambda_star = lam + alpha_star.* w;
                phi_star = Gauss_integral3D(lambda_star, 'f');
                fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi_star); 
            end
            return;
        end
        if phid_a >= 0
        
            [alpha_star, errCode] = zoom(lam, w, alpha(i), alpha(i-1), verbose);
            if verbose
                lambda_star = lam + alpha_star.* w;
                phi_star = Gauss_integral3D(lambda_star, 'f');
                fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi_star); 
            end
            return;
        end
        
        % choose next alpha
        alpha(i+1) = alpha(i) * multiplier;
        phi_last = phi_a;
    end
    
    % Exceeds alpha_max
    alpha_star = -9;
    if verbose
        lambda_star = lam + alpha_star.* w;
        phi_star = Gauss_integral3D(lambda_star, 'f');
        fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi_star); 
    end
    errCode = -1;
end


function [alpha_star, errCode] = zoom(lam, w, alpha_lo, alpha_hi, verbose)

    c1 = 0.0001;
    c2 = 0.9;
    phi_0 = Gauss_integral3D(lam, 'f');
    phid_0 = Gauss_integral3D(lam, 'g')' * w;

    while abs(alpha_hi - alpha_lo) > 0.00001
        
        % Cubic interpolation
        lambda_lo = lam + alpha_lo .* w;
        lambda_hi = lam + alpha_hi .* w;
        phi_lo = Gauss_integral3D(lambda_lo, 'f');
        phi_hi = Gauss_integral3D(lambda_hi, 'f');
        phid_lo = Gauss_integral3D(lambda_lo, 'g')' * w;
        phid_hi = Gauss_integral3D(lambda_hi, 'g')' * w;
        d1 = phid_lo + phid_hi - 3 * (phi_lo - phi_hi) / (alpha_lo - alpha_hi);
        d2 = sign(alpha_hi - alpha_lo) * sqrt(d1^2 - phid_lo * phid_hi);
        alpha = alpha_hi - (alpha_hi - alpha_lo) * (phid_hi + d2 - d1) / (phid_hi - phid_lo + 2*d2);
        % print
        if verbose, fprintf('** %f %f -> %f **\n', alpha_hi, alpha_lo, alpha);end
        % interior check
        if imag(alpha) ~= 0 || (alpha <= alpha_hi && alpha <= alpha_lo) || (alpha >= alpha_hi && alpha >= alpha_lo)
        
            if phi_lo <= phi_hi
                
                if alpha_lo==0
                    alpha_star = 0.0001;
                else 
                    alpha_star = alpha_lo;
                end
                fprintf('No interior solution -> %f\n', alpha_star);
                errCode = 0;
                return;
            else
                if alpha_hi==0
                    alpha_star = 0.0001;
                else 
                    alpha_star = alpha_hi;
                end
                fprintf('No interior solution -> %f\n', alpha_star);
                errCode = 0;
                return;
            end
        end
        
        % Zoom logic
        % evaluate phi(alpha)
        
        phi_a = Gauss_integral3D(lam + alpha.*w, 'f');
        if phi_a > (phi_0 + c1 * alpha * phid_0) || (phi_a >= phi_lo)

            alpha_hi = alpha;
        else

            % evaluate phid(alpha)
            phid_a = Gauss_integral3D(lam + alpha.*w, 'g')' * w;
            if abs(phid_a) <= -c2 * phid_0

                alpha_star = alpha;
                errCode = 0;
                return;
            end
            if phid_a * (alpha_hi - alpha_lo) >= 0

                alpha_hi = alpha_lo;
            end
            alpha_lo = alpha;
        end
    end

    % Exceed limitation
    alpha_star = 0.0001;
    errCode = 0;
end