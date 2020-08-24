%% Flexible dimension version. % Final version !
% Only need to edit this part.
clear; clc;


% source complexes of the reactions. (each complex could be repeated.)
K = 4; 
d = 2;

% d: the number of species
% K: the number of reactions.
syms n [1 d] integer
syms gg positive
ncell = sym2cell(n);
lambda(n) = sym(zeros(K,1));
lambda_cell = sym2cell(formula(lambda));

% initialization of the propensity functions
% complexes = [0 0;2 0;1 1;0 1];
% alpha = sym([3,1,2,0.5]);
% lambda_cell{1}(n) = alpha(1);
% lambda_cell{2}(n) = alpha(2) * n(1) * (n(1) - 1);
% lambda_cell{3}(n) = alpha(3) * n(1) * n(2);
% lambda_cell{4}(n) = alpha(4) * n(2);

% test set 2
% complexes = [0 0;1 0;0 1];
% alpha = sym([3,1,2]);
% syms alpha [1 K] positive
% lambda_cell{1}(n) = alpha(1) * (n(1) + gg);
% lambda_cell{2}(n) = alpha(2) * n(1) * (n(1) + gg);
% lambda_cell{3}(n) = alpha(3) * n(2) * (n(1) + gg);

% test set 3 - Fig. 2A
% alpha = sym([3,1,2]);
complexes = [0 0;1 0;0 1;0 1];
syms alpha [1 K] positive
lambda_cell{1}(n) = alpha(1);
lambda_cell{2}(n) = alpha(2) * n(1) * (n(1) - 1);
lambda_cell{3}(n) = alpha(3) * n(2);
lambda_cell{4}(n) = alpha(4) * n(2);

[a_list_tmp, b_list_tmp, elementary_basis] = CRN_find_elemtary_path(complexes); % the row vectors of H form a basis.
F = CRN_find_elementary_function(complexes, lambda_cell, alpha);
start_point = 3*ones(1,d);
elementary_coordinates = CRN_solve_sym_linear(elementary_basis, start_point);
coord_struct = struct2cell(elementary_coordinates);
theta(ncell{:}) = simplify(CRN_theta_construction(start_point, coord_struct, elementary_basis, F));


%% Compute complex balanced equilbrium.

syms ca cb k1 k2 k3 k4
eqn1 = k1 == cb * k3;
eqn2 = k1+k4*cb == ca * k2;
he = solve(eqn1,eqn2,ca,cb);

syms ca cb k1 k2 k3 k4
eqn1 = k1 == ca *cb* k4;
eqn2 = k1+k3*ca*cb == ca * k2;
sol1 = solve(eqn1,eqn2,ca,cb);
sol1.ca
sol1.cb
