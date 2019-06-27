% m = load('data/m1.txt');
m = [29.7352; 0.1439; 0.1257; 0.1117; 0.1004];

% the integrand part in 'f'
f_int = @(x, lam, n) exp(lam'*power(x, 0:n-1)');
% build the core function
f = @(lam, m, n) integral_impl(@(x) f_int(x, lam, n), 0, 1) - lam'*m;
% the jacobian function
g_int = @(i, j, x, lam, n) x.^(i+j).*exp(lam'*power(x, 0:n-1)');
g = @(i, j, lam, n) integral_impl(@(x) g_int(i, j, x, lam, n), 0, 1);
% partial derivative function
p_int = @(i, x, lam, n) x.^i.*exp(lam'*power(x, 0:n-1)');
p = @(i, lam, m, n) integral_impl(@(x) p_int(i, x, lam, n), 0, 1) - m(i+1);


% l0 = -ones(101,1);
l0 = [-1; -1; -1; -1; -1];
prec= 1e-6;
maxIter = 1000;
[minf, lmin] = BFGS(f, p, l0, prec, m, maxIter);
% [minf, lmin] = BFGS_cut(f, p, l0, prec, m, maxIter);
minf
