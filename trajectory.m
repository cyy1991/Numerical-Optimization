% Plot the trajectory

syms x_
m = [0.5, 0.1];

if (errCode ~= 0)
    disp("Warning: The algorithm has failed and the result is not credible. Comment me to override.");
    return;
end

maxX = max(xhist(:, 1));
minX = min(xhist(:, 1));
marX = 0.1*(maxX - minX);
maxY = max(xhist(:, 2));
minY = min(xhist(:, 2));
marY = 0.1*(maxY - minY);
[xx, yy] = meshgrid(minX-marX:(maxX-minX)/101:maxX+marX,...
                    minY-marY:(maxY-minY)/101:maxY+marY);
% zz = int(exp(xx + yy.*x_), x_, 0, 1) - (xx.*m(1) + yy.*m(2));
% zz = double(zz);
f = @(x, y) (exp(x+y)-exp(x))./y - m(1).*x - m(2).*y;
zz = zeros(size(xx));
for i = 1:size(xx, 1)
    for j = 1:size(xx, 2)
    
        zz(i, j) = f(xx(i, j), yy(i, j));
    end
end

maxZ = max(fhist);
minZ = min(fhist);
zz(zz > 1.1*maxZ-0.1*minZ) = 1.1*maxZ-0.1*minZ;

surf(xx, yy, zz)

% The following code is used to plot the trajectory of iterations
% and requires 'xhist' 'fhist' defined
hold on
%scatter3(xhist(1:70, 1), xhist(1:70, 2), fhist(1:70), 'red');
 for i = 1:length(fhist)
    text(xhist(i, 1), xhist(i, 2), fhist(i)+0.01*(maxZ-minZ), num2str(i), 'Color', 'r', 'FontSize', 12);
 end
