% Surf function wuth int

syms x_
m = [0.1, 5];

[xx, yy] = meshgrid(-20:1:3, -10:1:4);
zz = int(exp(xx + yy.*x_), x_, 0, 1) - (xx.*m(1) + yy.*m(2));
zz = double(zz);

surf(xx, yy, zz)
