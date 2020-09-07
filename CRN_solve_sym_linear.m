function solutions = CRN_solve_sym_linear(basis, start_point)
% basis: m * d matrix whose row vectors are basis of an augmented state space
% determined by Hermite normal form. i.e., zeta_i's
% 'augmented' means that we consider not only stoichiometric vector but
% also difference between any two source complexes.
% The descrepancy by this augmentation appears especially when there are
% multiple linkage classes.
% start_point: d*1 vector, start point that path from to an arbitrary state n.

% s: dimension of the augmented state space. 
% d: the number of species

Z = null(basis, 'r');

d = size(basis,2);
s = rank(basis);
basis = basis(1:s, :);
n_vector = sym('n', [1 d], 'integer');

assumeAlso((n_vector - start_point) * Z == 0) 
% this assumption forces n in the augmented state space.

c_vector = sym('c', [1 s], 'integer');
solutions = solve(c_vector * basis == n_vector - start_point, c_vector, 'ReturnConditions', true);

end
