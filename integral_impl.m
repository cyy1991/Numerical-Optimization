function [res] = integral_impl (f, start_, end_)
%% Different implementations of integral

    % res = default_integral(f, start_, end_);
    res = space_sample_integral(f, start_, end_, 1000);
    
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
end
