function cbe = CRN_compute_cbe(complexes, M_kappa)
% Compute complex balanced equilbrium (CBE) for the deterministic mass-action
% kinetic model with a given kinetic parameters kappa.

% data input
% complexes: d X |C| matrix whose column vectors are complex vector.
% note that unlike propensity factorization code, this 'complexes' variable
% does not contain repeated complexes. All the complexes have to be
% enumerated without overlapping.
% M_kappa: |C| X |C| matrix. (i,j)-entry represents the kinetic parameter of
% the reaction from the j-th complex to the i-th complex. A zero in the
% matrix represents there is no reaction. 
% M_kappa can be either numeric or symbolic matrix.

[d, num_C] = size(complexes);
syms c [d 1] positive
intensities = (prod(c.^complexes, 1))';
% intensities: num_C X 1 column vector representing the intensity functions of the
% reactions from each of complexes without the kinetic parameters.
% The underlying kinetics is assumed to be mass-action.

A_k = M_kappa;
for j = 1:num_C
    A_k(j, j) = -sum(M_kappa(:, j));
end
% Laplace matrix.

eqn1 = A_k * intensities == 0;
cbe = solve(eqn1,c);


%%% The below code is not necessary. The code is built for incorporating
%%% the characteristics of CBE. (solutions of log-linear systems)

reactionTF = zeros(num_C, num_C);

for i = 1:num_C
    for j = 1:num_C
        if ~isequal(M_kappa(i,j), sym(0))
            reactionTF(i,j) = 1;
        end
    end
end
num_R = sum(sum(reactionTF));

% stoichiometric matrix 
stoi_M = zeros(d, num_R);
k = 0;
for j = 1:num_C
    for i = 1:num_C
        if reactionTF(i,j) == 1
            k = k+1;
            stoi_M(:,k) = complexes(:,i) - complexes(:,j);
        end
    end
end

syms c0
conservation_law = null(stoi_M', 'r')' * c == c0;

end