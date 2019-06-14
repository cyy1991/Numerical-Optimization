% Plot the trajectory

syms x_
m = [0.5, 0.1];

[xx, yy] = meshgrid(min(xhist(:, 1))-1.05:0.1:max(xhist(:, 1))+1.05,...
                    min(xhist(:, 2))-1.05:0.1:max(xhist(:, 2))+1.05);
% zz = int(exp(xx + yy.*x_), x_, 0, 1) - (xx.*m(1) + yy.*m(2));
% zz = double(zz);
f = @(x, y) (exp(x+y)-exp(x))./y - m(1).*x - m(2).*y;
zz = zeros(size(xx));
for i = 1:size(xx, 1)
    for j = 1:size(xx, 2)
    
        zz(i, j) = f(xx(i, j), yy(i, j));
    end
end

zz(zz > max(fhist)*2) = max(fhist)*2;

surf(xx, yy, zz)

% The following code is used to plot the trajectory of iterations
% and requires 'xhist' 'fhist' defined
hold on
%scatter3(xhist(1:70, 1), xhist(1:70, 2), fhist(1:70), 'red');
 for i = 1:length(fhist)
    text(xhist(i, 1), xhist(i, 2), fhist(i)+0.1, num2str(i), 'Color', 'r', 'FontSize', 12);
 end
