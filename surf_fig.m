% Surf function wuth int

syms x_
m = [0.5, 0.1];

[xx, yy] = meshgrid(-3:0.3:3, -7:0.3:2);
% zz = int(exp(xx + yy.*x_), x_, 0, 1) - (xx.*m(1) + yy.*m(2));
% zz = double(zz);
f = @(x, y) (exp(x+y)-exp(x))./y - m(1).*x - m(2).*y;
zz = zeros(size(xx));
for i = 1:size(xx, 1)
    for j = 1:size(xx, 2)
    
        zz(i, j) = f(xx(i, j), yy(i, j));
    end
end


zz(zz > 10) = 10;

surf(xx, yy, zz)
