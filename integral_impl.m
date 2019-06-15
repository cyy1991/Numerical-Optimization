function [res] = integral_impl (f, start_, end_)
%% Different implementations of integral

    % res = default_integral(f, start_, end_);
    % res = space_sample_integral(f, start_, end_, 1000);
    res = romberg_integral(f, start_, end_, 10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Implementations start here
    function [res] = default_integral (f, start_, end_)
        
        % Unknown bug when integrating f
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
            
            sum = 0;
            for k = 1:2^(i-2)
                sum = sum + f(start_ + (k-0.5)*h);
            end
            R(2, 1) = 1/2 * (R(1, 1) + h*sum);

            for j = 2:i
                R(2, j) = R(2, j-1) + (R(2, j-1)-R(1, j-1)) / (4^(j-1)-1);
            end

            h = h/2;
            R(1, :) = R(2, :);
        end

        res = R(2, n);
    end
end
