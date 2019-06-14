function [res] = integral_impl (f, start_, end_)

    % res = default_integral(f, start_, end_);
    res = space_sample_integral(f, start_, end_, 200);
    
    %% Implementations
    function [res] = default_integral (f, start_, end_)
        
        % Unknown bug when integral f
        res = integral(f, start_, end_);
    end

    function [res] = space_sample_integral (f, start_, end_, n)
    % n is the number of evaluations
        
        xx = linspace(start_, end_, n);
        yy = arrayfun(f, xx);
        res = sum(yy) / n;
    end
end
