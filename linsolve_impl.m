function [x] = linsolve_impl (A, b)
%% Solve linear system Ax=b

    x = default_linsolve(A, b);
    % x = guass_scaled_partial(A, b);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Implementations start here
    function [x] = default_linsolve(A, b)
    
        x = linsolve(A, b);
    end

    function [x] = guass_scaled_partial(A, b)
    
        n = length(b);
        v = zeros(n, 1);
    end
end
