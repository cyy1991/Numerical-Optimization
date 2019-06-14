% Surf function wuth int

syms x_
m = [0.5, 0.1];

[xx, yy] = meshgrid(-3:0.5:3, -7:0.5:2);
zz = int(exp(xx + yy.*x_), x_, 0, 1) - (xx.*m(1) + yy.*m(2));
zz = double(zz);
% zz_simp = (exp(xx+yy)-exp(xx))./yy - 0.5.*xx - yy;
zz(zz > 10) = 10;

surf(xx, yy, zz)

% The following code is used to plot the trajectory of iterations
% and requires 'xhist' 'fhist' defined
hold on
%scatter3(xhist(1:70, 1), xhist(1:70, 2), fhist(1:70), 'red');
 for i = 1:7
    text(xhist(i, 1), xhist(i, 2), fhist(i)+1, num2str(i), 'Color', 'r', 'FontSize', 12);
 end
