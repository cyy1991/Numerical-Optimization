% test for linsolve_impl

N = 20;
n = 1000;

err = zeros(n, 1);

timer_(-1);

for i = 1:N

    A = rand(n);
    y = rand([n, 1]);
    
    timer_(0);
    y1 = linsolve(A, y);
    timer_(1);
    y2 = linsolve_impl(A, y);
    timer_(2);
    err(i) = sum(y1 - y2);
end

sum(err)
timer_(-3);
