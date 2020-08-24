function [theta, n] = CRN_construct_theta(nmlzd_lambdas, complexes, path, b)
% nmlzd_lambdas: (R * 1) matrix contaning propensities / kappa for all reactions.
% path: (2*m) matrix whose column [aj, bj] represents a unit path, (nu_aj - nu_bj). 
% complexes: d * |C| matrix whose columns are complex vectors. 
% b: (d * 1) vector. Reference (or start) point of the path from b to n.

m = size(path,2);
theta_tmp = 1;
nj = b;
for jj = 1:m
    numerator_fun = nmlzd_lambdas{path(1,jj)};
    denominator_fun = nmlzd_lambdas{path(2,jj)};
    theta_tmp = numerator_fun(nj + complexes(:,path(1,jj)) - complexes(:,path(2,jj))) / ...
        denominator_fun(nj) * ...
        theta_tmp;
    nj = nj + complexes(:,path(1,jj)) - complexes(:,path(2,jj));
end
theta = theta_tmp;
n = nj;
end

% example setting;
nmlzd_lambdas = {@(x)1,@(x) x(1),@(x) x(2)};
complexes = [0,0;1,0;0,1]';
path = [2,1;2,1;2,1;2,1;2,1;2,1;2,1;2,1;3,2;3,2;3,2;3,2;3,2]';
b = [0;0];

syms n
p = sym(1);

for i = sym(1):n+3
    p = p * 1/i;
end

syms n i
symprod(i, i, 1,n)


M = [1 0; 1 1; 2 0; -2 1; -1 1; -1 0]



