function f_int = Gauss_integral3D(i1_i2_i3_lambda)
% INPUT: [i1, i2, i3, lambda] is a (d x 4) matrix
% directly implement the integral from 0 to 1, given lambda
% use five point Gauss in 3D, in total 125 evaluation point
% reach accuracy of at least 5-6 decimal point

% example: 
% Gauss_integral3D([1,0, 0, 0.2])
% Gauss_integral3D([3,3,0,0.01; 2,2,1,0.4;1,3,1,0.006])

% base function change of variable
each_term = @(m, t, s, i1, i2, i3, lambda) lambda * m.^i1 .* ((1-m).*t).^i2 .* ((1-m-(1-m).*t).*s).^i3; 

% generate zero_node/weights in 3D space
zero_nodes_base = 0.5* [-0.90618, -0.538469, 0, 0.538469, 0.90618 ] + 0.5;
w = [0.236927; 0.478629; 0.568889; 0.478629; 0.236927];
zero_nodes = zeros(125, 3);
weights = zeros(1,125);
for i=1:5
    for j=1:5
        for k=1:5
            zero_nodes(25*(i-1)+5*(j-1)+(k),: ) = zero_nodes_base([ i, j, k]);
            weights(25*(i-1)+5*(j-1)+(k) )  = w(i)*w(j)*w(k);
        end
    end
end
m = zero_nodes(:,1);
t = zero_nodes(:,2);
s = zero_nodes(:,3);
change_variable_constant = (1-m).*(1-m-(1-m).*t);
f_val = zeros(125,1);
[row_n, col_n] = size(i1_i2_i3_lambda);
for i = 1:row_n
    f_val = f_val + each_term( m, t, s,  i1_i2_i3_lambda(i,1), i1_i2_i3_lambda(i,2), i1_i2_i3_lambda(i,3),i1_i2_i3_lambda(i,4) );
end
f_val = exp(f_val) .* change_variable_constant ;
f_int = weights * f_val  / 8; 

end