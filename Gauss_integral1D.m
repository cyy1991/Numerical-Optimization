function f_int = Gauss_integral1D(lambda)
% directly implement the integral from 0 to 1, given lambda
% use three point Gauss 
x = 0.5* [-0.90618, -0.538469, 0, 0.538469, 0.90618 ] + 0.5;
w = [0.236927; 0.478629; 0.568889; 0.478629; 0.236927];
n = length(lambda);
f = [];
for i=1:5
    f = [f, exp(lambda*power(x(i), 0:n-1)')];
end
f_int = 0.5* f * w;
end