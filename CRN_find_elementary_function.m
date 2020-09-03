function F = CRN_find_elementary_function(complexes, lambda_cell, kappa)
% complexes: d * K matrix whose column vectors are the source complexs of
% reactions.
[d, ~] = size(complexes);

syms n [d 1] integer
ncell = sym2cell(n);

[a_list_tmp, b_list_tmp, H] = CRN_find_elemtary_path(complexes);
% the row vectors of H form a basis.

number_of_elementary_path = rank(H);
length_of_elementary_path = sum(a_list_tmp ~= 0, 2); % compute row sum.

F_tmp(n) = sym(zeros(number_of_elementary_path, 1));
F = sym2cell(formula(F_tmp));

for l = 1:number_of_elementary_path
    a_list = a_list_tmp(l, 1:length_of_elementary_path(l));
    b_list = b_list_tmp(l, 1:length_of_elementary_path(l));
    
    trace_of_current_path = cumsum(complexes(:,a_list) - complexes(:,b_list), 2);

    for j = 1:length_of_elementary_path(l)
        tmp_arg_n_a = sym2cell(cellfun(@sum, ncell) + trace_of_current_path(:,j));
        
        if j == 1
            F{l}(ncell{:}) = (lambda_cell{a_list(j)}(tmp_arg_n_a{:})/kappa(a_list(j))) / ...
                (lambda_cell{b_list(j)}(ncell{:}) / kappa(b_list(j)));
        else
            tmp_arg_n_b = sym2cell(cellfun(@sum, ncell) + trace_of_current_path(:,j-1));
            F{l}(ncell{:}) = F{l}(ncell{:}) * (lambda_cell{a_list(j)}(tmp_arg_n_a{:})/kappa(a_list(j))) / ...
                (lambda_cell{b_list(j)}(tmp_arg_n_b{:}) / kappa(b_list(j)));
        end
    end
end

end