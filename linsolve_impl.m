function [x] = linsolve_impl (A, b, choice)
%% Solve linear system Ax=b

    if (nargin < 3)
    
        choice = 2;
    end
    if (choice == 1)
        x = default_linsolve(A, b);
    else
        x = guass_scaled_partial(A, b);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Implementations start here
    function [x] = default_linsolve(A, b)
    
        x = linsolve(A, b);
    end

    function [x] = guass_scaled_partial(A, b)

        % Guass with partial pivoting
        [n, m] = size(b);
        v = 1 : n;  % Row swapping recorder
        s = zeros(n, 1);

        for i = 1:n

            s(i) = max(abs(A(i, :)));
        end

        for k = 1:n-1

            rmax = 0;
            for i = k:n

                r = abs(A(v(i), k) / s(v(i)));
                if r > rmax
                    rmax = r;
                    t = i;
                end
            end

            temp = v(k);
            v(k) = v(t);
            v(t) = temp;

            for i = k+1:n

                mul = A(v(i), k) / A(v(k), k);
                A(v(i), k) = mul;  % Only serve the purpose of recorder
                for j = k+1:n
                    A(v(i), j) = A(v(i), j) - mul*A(v(k), j);
                end
            end
        end

        x = zeros(n, m);

        % Back substitution
        for k = 1:m

            for j = 1:n-1
                for i = j+1:n
                    b(v(i), k) = b(v(i), k) - A(v(i), j)*b(v(j), k);
                end
            end

            x(n, k) = b(v(n), k) / A(v(n), n);

            for i = n-1:-1:1

                z = b(v(i), k);
                for j = i+1:n
                    z = z - A(v(i), j)*x(j, k);
                end
                x(i, k) = z / A(v(i), i);
            end
        end

    end
end
