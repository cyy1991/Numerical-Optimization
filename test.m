% test for linsolve

N = 1;
n = 4;

for i = 1:N

    A = rand(n);
    y = rand([n, 1]);
    
    y1 = linsolve(A, y)
    y2 = linsolve_impl(A, y)
    sum(y1 - y2)
end
