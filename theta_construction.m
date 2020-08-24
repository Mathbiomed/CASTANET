%% Flexible dimension version. % Final version !
% Only need to edit this part.
clear; clc;
d = 2;
syms n [1 d] integer
K = 4; % the number of reactions.
ncell = sym2cell(n);

% initialization of the propensity functions
alpha = sym([3,1,2,0.5]);

lambda(n) = sym(zeros(K,1));
lambda_cell = sym2cell(formula(lambda));

lambda_cell{1}(n) = alpha(1);
lambda_cell{2}(n) = alpha(2) * n(1) * (n(1) - 1);
lambda_cell{3}(n) = alpha(3) * n(1) * n(2);
lambda_cell{4}(n) = alpha(4) * n(2);

complexes = [0 0;2 0;1 1;0 1];
% source complexes of the K reactions. (a complex could be repeated.)

[a_list_tmp, b_list_tmp, H] = CRN_find_elemtary_path(complexes);
% the row vectors of H form a basis.

number_of_elementary_path = rank(H);
length_of_elementary_path = sum(a_list_tmp ~= 0, 2); % compute row sum.

F_tmp(n) = sym(zeros(number_of_elementary_path, 1));
F = sym2cell(formula(F_tmp));

for l = 1:number_of_elementary_path
    a_list = a_list_tmp(l, 1:length_of_elementary_path(l));
    b_list = b_list_tmp(l, 1:length_of_elementary_path(l));
    
    trace_of_current_path = cumsum(complexes(a_list,:) - complexes(b_list,:));

    for j = 1:length_of_elementary_path(l)
        tmp_arg_n_a = sym2cell(cellfun(@sum, ncell) + trace_of_current_path(j,:));
        
        if j == 1
            F{l}(ncell{:}) = (lambda_cell{a_list(j)}(tmp_arg_n_a{:})/alpha(a_list(j))) / ...
                (lambda_cell{b_list(j)}(n1,n2) / alpha(b_list(j)));
        else
            tmp_arg_n_b = sym2cell(cellfun(@sum, ncell) + trace_of_current_path(j-1,:));
            F{l}(ncell{:}) = F{l}(ncell{:}) * (lambda_cell{a_list(j)}(tmp_arg_n_a{:})/alpha(a_list(j))) / ...
                (lambda_cell{b_list(j)}(tmp_arg_n_b{:}) / alpha(b_list(j)));
        end
    end
end

elementary_basis = H(1:rank(H), :);
start_point = zeros(1,d);
elementary_coordinates = CRN_solve_sym_linear(elementary_basis, start_point);
coord_struct = struct2cell(elementary_coordinates);

stopover = sym(start_point);

theta = sym(1);
syms j

for l = 1:number_of_elementary_path
    tmp = stopover + (j-1) * elementary_basis(l,:);
    stopover_cell = num2cell(tmp);
    theta = theta * symprod(F{l}(stopover_cell{:}), j, 1, coord_struct{l});
    if l > 1
        stopover = stopover + coord_struct{l-1} * elementary_basis(l-1,:);
    end
end



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
