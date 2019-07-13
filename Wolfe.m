function [res] = Wolfe(f, g, x, p, alpha, c1, c2)
% Wolfe Conditions

    if (nargin <= 5)
    
        c1 = 0.0001;
        c2 = 0.9;
    end

    f_x = f(x);
    f_xp = f(x + alpha .* p);
    g_x = g(x);
    g_xp = g(x + alpha .* p);
    
    res = (f_xp <= f_x + c1 .* alpha .* p' * g_x) & ...
          (- p' * g_xp <= - c2 .* p' * g_x);
end
