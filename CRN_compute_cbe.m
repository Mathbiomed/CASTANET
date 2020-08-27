function cbe = CRN_compute_cbe(complexes, M_kappa)
% Compute complex balanced equilbrium (CBE) for the deterministic mass-action
% kinetic model with a given kinetic parameters kappa.

% data input
% complexes: |C| X d matrix whose row vectors are complex vector.
% note that unlike propensity factorization code, this 'complexes' variable
% does not contain repeated complexes. All the complexes have to be
% enumerated without overlapping.
% M_kappa: |C| X |C| matrix. (i,j)-entry represents the kinetic parameter of
% the reaction from the j-th complex to the i-th complex. A zero in the
% matrix represents there is no reaction. 
% M_kappa can be either numeric or symbolic matrix.


% test case
% syms kappa1 kappa2 kappa3 kappa4 positive
% complexes = [1 1 0; 0 2 0; 0 0 1];
% M_kappa = [0, kappa2, kappa4; kappa1, 0, 0; 0, kappa3, 0];

[num_C, d] = size(complexes);
sources = prod(c.^complexes, 2);
% sources: num_C X 1 column vector representing the intensity functions of the
% reactions from each of complexes without the kinetic parameters.
syms c [1 d] positive

A_k = M_kappa;
for j = 1:num_C
    A_k(j, j) = -sum(M_kappa(:, j));
end

eqn1 = A_k * sources == 0;
sol1 = solve(eqn1,c);



reactionTF = zeros(num_C, num_C);

for i = 1:num_C
    for j = 1:num_C
        if ~isequal(M_kappa(i,j), sym(0))
            reactionTF(i,j) = 1;
        end
    end
end
num_R = sum(sum(reactionTF));
% A_k is the matrix satisfying A_k * sources = 0 implies c is a CBE.
% stoichiometric matrix 
stoi_M = zeros(num_R, d);
k = 0;
for j = 1:num_C
    for i = 1:num_C
        if reactionTF(i,j) == 1
            k = k+1;
            stoi_M(k,:) = complexes(i,:) - complexes(j,:);
        end
    end
end

syms c0

conservation_law = null(stoi_M, 'r')' * c' == c0;


% examples - this part should be removed after complete the code.
syms ca cb k1 k2 k3 k4
eqn1 = k1 == cb * k3;
eqn2 = k1+k4*cb == ca * k2;
sol1 = solve(eqn1,eqn2,ca,cb);

syms ca cb k1 k2 k3 k4
eqn1 = k1 == ca *cb* k4;
eqn2 = k1+k3*ca*cb == ca * k2;
sol1 = solve(eqn1,eqn2,ca,cb);
sol1.ca
sol1.cb

end