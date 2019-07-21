function [res] = integral_impl (f, start_, end_, mandatory)
%% Different implementations of integral

    persistent prec count

    % res = default_integral(f, start_, end_);
    % res = space_sample_integral(f, start_, end_, 1000);
    % res = romberg_integral(f, start_, end_, 20);

    if nargin > 3

        if mandatory == 0

            prec = 4;
            count = 0;
        end
        if mandatory == 5
        
            prec = 5;
        end
        return;
    end

    count = count + 1;
    if prec == 1

        res = guassian_quadrature_100(f);
        if mod(count, 1000) == 0
        
            count = 1;
            fprintf('Cross validating integral accuracy...\n');
            res_auth = guassian_quadrature_250(f);
            if abs(res-res_auth) > 1
        
                prec = prec + 1;
                fprintf('-------------------------\n|Integral precision = %d\n-------------------------\n', prec);
                res = res_auth;
                return;
            end
        end
    elseif prec == 2
        
        res = guassian_quadrature_250(f);
        if mod(count, 100) == 0
        
            count = 1;
            fprintf('Cross validating integral accuracy...\n');
            res_auth = guassian_quadrature_500(f);
            if abs(res-res_auth) > 1
        
                prec = prec + 1;
                fprintf('-------------------------\n|Integral precision = %d\n-------------------------\n', prec);
                res = res_auth;
                return;
            end
        end
    elseif prec == 3
        
        res = guassian_quadrature_500(f);
        if mod(count, 50) == 0
        
            count = 1;
            fprintf('Cross validating integral accuracy...\n');
            res_auth = guassian_quadrature_1000(f);
            if abs(res-res_auth) > 1
        
                prec = prec + 1;
                fprintf('-------------------------\n|Integral precision = %d\n-------------------------\n', prec);
                res = res_auth;
                return;
            end
        end
    elseif prec == 4

        res = guassian_quadrature_1000(f);
    elseif prec == 5

        res = guassian_quadrature_10000(f);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Implementations start here
    function [res] = default_integral (f, start_, end_)

        % Unknown bug when integrate f
        res = integral(f, start_, end_);
    end

    function [res] = space_sample_integral (f, start_, end_, n)
    % n is the number of evaluations

        xx = linspace(start_, end_, n);
        yy = arrayfun(f, xx);
        res = sum(yy) / n;
    end

    function [res] = romberg_integral (f, start_, end_, n)
    % n is the size of Romberg table

        % initialization
        h = end_ - start_;
        R = zeros(2, n);
        R(1, 1) = h/2 * (f(start_) + f(end_));

        for i = 2:n
            
            sum_ = sum(arrayfun(f, start_ + ((1:2^(i-2))-0.5)*h));
            R(2, 1) = 1/2 * (R(1, 1) + h*sum_);

            for j = 2:i
                R(2, j) = R(2, j-1) + (R(2, j-1)-R(1, j-1)) / (4^(j-1)-1);
            end

            h = h/2;
            R(1, :) = R(2, :);
        end

        res = R(2, n);
    end

    function [res] = adaptive_romberg_integral(f, start_, end_, n)
    
        h = end_ - start_;
        R = zeros(2, n);
        R(1, 1) = h/2 * (f(start_) + f(end_));

        record = zeros(n, 1);
        for i = 2:n
            
            sum_ = sum(arrayfun(f, start_ + ((1:2^(i-2))-0.5)*h));
            R(2, 1) = 1/2 * (R(1, 1) + h*sum_);

            for j = 2:i
                R(2, j) = R(2, j-1) + (R(2, j-1)-R(1, j-1)) / (4^(j-1)-1);
            end

            h = h/2;
            R(1, :) = R(2, :);
            record(i) = R(2, i);
            if i > 9 && abs(record(i) - record(i-1)) < 1
            
                res = record(i);
                % fprintf('[%d/%d]', i, n);
                return;
            end
        end

        res = R(2, n);
    end

    function [res] = guassian_quadrature_100(f)
    % n is the number of nodes
    
        global w_100 x_100

        ff = arrayfun(f, x_100);
        res = 0.5 .* ff' * w_100;
    end

    function [res] = guassian_quadrature_250(f)
    % n is the number of nodes
    
        global w_250 x_250
        
        ff = arrayfun(f, x_250);
        res = 0.5 .* ff' * w_250;
    end

    function [res] = guassian_quadrature_500(f)
    % n is the number of nodes
    
        global w_500 x_500
        
        ff = arrayfun(f, x_500);
        res = 0.5 .* ff' * w_500;
    end

    function [res] = guassian_quadrature_1000(f)
    % n is the number of nodes
        
        global w_1000 x_1000
        
        ff = arrayfun(f, x_1000);
        res = 0.5.* ff' * w_1000;
    end

    function [res] = guassian_quadrature_10000(f)
    % n is the number of nodes
        
        global w_10000 x_10000
        
        ff = arrayfun(f, x_10000);
        res = 0.5.* ff' * w_10000;
    end
end
