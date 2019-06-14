% Surf function wuth int

syms x_
m = [0.5, 1];

[xx, yy] = meshgrid(-5:0.5:5, -5:0.5:6);
zz = int(exp(xx.*x_ + yy.*x_.^2), x_, 0, 1) - (xx.*m(1) + yy.*m(2));
zz = double(zz);
zz(zz > 30) = 30;

surf(xx, yy, zz)
