function factorization_TF = CRN_check_factorization_condition(complexes, lambda_cell, theta, kappa, ncell)
% Check the factorization conditions for all the reactions.
% compute the function omega using the condition of the first reaction.
[~, K] = size(complexes);
% syms n [d 1] integer
% ncell = sym2cell(n);
tmp_arg_n_nu1 = sym2cell(cellfun(@sum, ncell) + complexes(:,1));
omega(ncell{:}) = lambda_cell{1}(tmp_arg_n_nu1{:}) / (kappa(1) *theta(tmp_arg_n_nu1{:})); 

% check all the condtions.
factorization_TF = zeros(1,K);
for k = 1:K
    tmp_arg_n_nuk = sym2cell(cellfun(@sum, ncell) - complexes(:,k));
    factorization_TF(k) = isAlways(lambda_cell{k} == kappa(k) * theta(ncell{:}) * omega(tmp_arg_n_nuk{:}));
end

% if prod(factorization_TF) == 1
%     disp("All the factorization conditions hold!");
% else
%     disp("Some of factorization conditions do not hold!");
%     disp("Please check the output to find the reactions that do not satisfy the conditions");
% end

end