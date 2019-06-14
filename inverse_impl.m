function [A_inv] = inverse_impl (A)
%% Different implementations of inverse matrix

    A_inv = default_inverse(A);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Implementations start here
    function [A_inv] = default_inverse(A)
    
        A_inv = inv(A);
    end
end
