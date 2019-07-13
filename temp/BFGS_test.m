% m = load('data/m1.txt');
m = [29.7352; 0.1439; 0.1257; 0.1117; 0.1004];



% l0 = -ones(101,1);
l0 = [-1; -1; -1; -1; -1];
prec= 1e-6;
maxIter = 1000;
[minf, lmin] = BFGS(f, p, l0, prec, m, maxIter);
% [minf, lmin] = BFGS_cut(f, p, l0, prec, m, maxIter);
minf
