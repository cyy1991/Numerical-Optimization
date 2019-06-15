function [x] = linsolve_impl (A, b)
%% Solve linear system Ax=b

    % x = default_linsolve(A, b);
    x = guass_scaled_partial(A, b);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Implementations start here
    function [x] = default_linsolve(A, b)
    
        x = linsolve(A, b);
    end

   function [x] = guass_scaled_partial(A, b)
        
        % guass with partial pivoting
        [n, m] = size(b);
        v = zeros(n, 1);
        s = zeros(n, 1);
        
        for i = 1:n
            v(i) = i;
            s(i) = max(abs(A(i,:)));
        end
        
        for k = 1:n-1
            rmax = 0;
            for i  = k:n
                r = abs(A(v(i),k)/s(v(i)));
                if r > rmax
                    rmax = r;
                    t = i;
                end
            end
            
            temp = v(k);
            v(k) = v(t);
            v(t) = temp;
            
            for i = k+1:n
                mul = A(v(i),k)/A(v(k),k);
                A(v(i),k) = mul;
                for j = k+1:n
                    A(v(i),j) = A(v(i),j)-mul*A(v(k),j);
                end
            end
        end
        
        x = zeros(n,m);
        
        % back substitution
        for k = 1:m
            for j = 1:n-1
                for i = j+1:n
                    b(v(i),k) = b(v(i),k) - A(v(i),j)*b(v(j),k);
                end
            end
            x(n,k) = b(v(n),k)/A(v(n),n);
            for i = n-1:-1:1
                z = b(v(i),k);
                for j = i+1:n
                    z = z-A(v(i),j)*x(j,k);
                end
                x(v(i),k) = z/A(v(i),i);
            end
        end
                
    end
end
