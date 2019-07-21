function [alpha_star, errCode] = LineSearch(phi, phid, w, it, verbose)

    if nargin<3, verbose = false;end
    c1 = 0.0001;
    c2 = 0.9;

    maxIt = 101;
    multiplier = 1.1482;
    alpha = zeros(maxIt+1, 1);
    alpha(1) = 0;
    alpha(2) = 1;
    % alpha_max = multiplier^(maxIt-1)
    errCode = 0;

    phi_0 = phi(0);
    phid_0 = phid(0);
    phi_last = 0;

    %if verbose, fprintf('------ iter: %d ------\n', it);end
    temp = phi(alpha(2));
    while isnan(temp) || (temp > 10e100)
        
        alpha(2) = alpha(2) / 10;
        temp = phi(alpha(2));
        if alpha(2) < 10e-30
            alpha(2) = 1;
            break;
        end
        if verbose, fprintf('{>>}');end
    end
    while isnan(temp) || (temp > 10e100)
        
        alpha(2) = alpha(2) + 0.05;
        temp = phi(alpha(2));
        if verbose, fprintf('{<<}');end
    end
    for i = 2:maxIt
    
        % evaluate phi(alpha)
        phi_a = phi(alpha(i));
        if verbose, fprintf('alpha:%f  phi:%f\n', alpha(i), phi_a);end
        if (phi_a > phi_0 + c1 * alpha(i) * phid_0) || (i > 2 && phi_a >= phi_last)

            [alpha_star, errCode] = zoom(phi, phid, w, alpha(i-1), alpha(i), verbose);
            if verbose, fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi(alpha_star)); end
            return;
        end
        % evaluate phid(alpha)
        phid_a = phid(alpha(i));
        if abs(phid_a) <= -c2 * phid_0
        
            alpha_star = alpha(i);
            if verbose, fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi(alpha_star)); end
            return;
        end
        if phid_a >= 0
        
            [alpha_star, errCode] = zoom(phi, phid, w, alpha(i), alpha(i-1), verbose);
            if verbose, fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi(alpha_star)); end
            return;
        end
        
        % choose next alpha
        alpha(i+1) = alpha(i) * multiplier;
        phi_last = phi_a;
    end
    
    % Exceeds alpha_max
    alpha_star = -9;
    if verbose, fprintf('Wolfe Alpha: %8.10f, Phi: %f\n', alpha_star,phi(alpha_star)); end
    errCode = -1;
end


function [alpha_star, errCode] = zoom(phi, phid, w, alpha_lo, alpha_hi, verbose)

    c1 = 0.0001;
    c2 = 0.9;
    phi_0 = phi(0);
    phid_0 = phid(0);

    while abs(alpha_hi - alpha_lo) > 10e-10

        % Cubic interpolation
        phi_lo = phi(alpha_lo);
        phi_hi = phi(alpha_hi);
        phid_lo = phid(alpha_lo);
        phid_hi = phid(alpha_hi);
        d1 = phid_lo + phid_hi - 3 * (phi_lo - phi_hi) / (alpha_lo - alpha_hi);
        d2 = sign(alpha_hi - alpha_lo) * sqrt(d1^2 - phid_lo * phid_hi);
        alpha = alpha_hi - (alpha_hi - alpha_lo) * (phid_hi + d2 - d1) / (phid_hi - phid_lo + 2*d2);
        % print
        if verbose, fprintf('** %f %f -> %f **\n', alpha_lo, alpha_hi, alpha);end
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
        phi_a = phi(alpha);
        if phi_a > (phi_0 + c1 * alpha * phid_0) || (phi_a >= phi_lo)

            alpha_hi = alpha;
        else

            % evaluate phid(alpha)
            phid_a = phid(alpha);
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
    alpha_star = 1000 / norm(w);
    fprintf('[ Warning: LineSearch stage II exceeds limitation. ]\n');
    errCode = -2;
end
