function [a_list, b_list, H] = CRN_find_elemtary_path(complexes)
% complexes: K * d matrix whose row vectors are the source complexs of
% reactions. Note that two of row vectors could be the same in case that
% the two reactions have a common source complex.
% K is the number of reactions, d is the number of species.
% Output H: the row vectors of H forms a basis with elementary paths.

[K, d] = size(complexes);

% The entries of A matrix is the difference between two source complexes.
A = zeros(K-1, d);
for k = 1:(K-1)
    A(k, :) = complexes(k+1,:) - complexes(1,:);
end
[U,H] = hermiteForm(A);

s = rank(A);
max_length_of_path = max(sum(abs(U),2));

a_list = zeros(s, max_length_of_path);
b_list = zeros(s, max_length_of_path);
% nonzero entries of l-th row vector of a_list is a_i indices from n to n + zeta_l;
% nonzero entries of l-th row vector of a_list is b_i indices from n to n + zeta_l;

for l = 1:s
    idx = 0;
    for k = 1:(K-1)
        for jj = 1:abs(U(l,k))
            idx = idx + 1;
            if U(l,k) > 0
                a_list(l, idx) = k+1;
                b_list(l, idx) = 1;
            elseif U(l,k) < 0
                a_list(l, idx) = 1;
                b_list(l, idx) = k+1;
            end
        end
    end
end
                
end
