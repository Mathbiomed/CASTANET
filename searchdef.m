% This is a main funtion for Network translation.

%% Initialization of all input variables and parameters
clear; clc;

% Fig. 1 example
sources = [0 0; 1 0; 0 1; 2 0]'; 
products = [1 0; 0 1; 0 0; 1 1]'; 
% (number of reaction) * (number of species) matrix containing the source and product complex vectors of reactions, respectively.
% 
% Fig. 2a example
% sources = [0 0; 2 0; 0 1; 0 1]'; 
% products = [1 0; 1 1; 0 0; 1 0]'; 

% % Fig. 2d example
% sources = [2 0 0; 0 1 0; 0 1 0; 0 0 1]'; 
% products = [1 1 0; 0 0 1;1 0 0;0 1 0]'; 

% % Fig. 2g example
% sources = [1 0; 1 1; 0 1]'; 
% products = [0 1; 0 2; 1 0]'; 
% 
% New example - necessary theta-omega
% sources = [0 0; 1 0; 0 1; 1 0; 2 0; 1 1]'; 
% products = [1 0; 0 1; 0 0; 2 0; 1 1 ; 1 0]'; 


[d, K] = size(sources);
%propensities = [];
% source complexes of the reactions. (each complex could be repeated.)
% K and d could be determined by [K, d] = size(complexes) after determining
% complexes but for the simplicity user would enter them here, manually.

% d: the number of species
% K: the number of reactions.
syms n [d 1] integer
% syms gg positive
ncell = sym2cell(n);
lambda(n) = sym(zeros(K,1));
lambda_cell = sym2cell(formula(lambda));
% tt = 1;
syms alpha [K 1] positive
% % Fig. 1 example

% lambda_cell{1}(n) = alpha(1);
% lambda_cell{2}(n) = alpha(2) * n(1);
% lambda_cell{3}(n) = alpha(3) * n(2);
% lambda_cell{4}(n) = alpha(4) * n(1) * (n(1) -1);

% 
% Fig. 2a example
% lambda_cell{1}(n) = alpha(1);
% lambda_cell{2}(n) = alpha(2) * n(1) * (n(1) -1);
% lambda_cell{3}(n) = alpha(3) * n(2);
% lambda_cell{4}(n) = alpha(4) * n(2);

% % Fig. 2d example
% lambda_cell{1}(n) = alpha(1) * n(1) * (n(1) -1);
% lambda_cell{2}(n) = alpha(2) * n(2);
% lambda_cell{3}(n) = alpha(3) * n(2);
% lambda_cell{4}(n) = alpha(4) * n(3);
% % 
% Fig. 2g example
% lambda_cell{1}(n) = alpha(1) * n(1);
% lambda_cell{2}(n) = alpha(2) * n(1) *n(2);
% lambda_cell{3}(n) = alpha(3) * n(2);
% assumeAlso(n(1)+n(2) == 20);
% % New example - necessary theta-omega
% syms g_const positive
% lambda_cell{1}(n) = alpha(1) * g_const;
% lambda_cell{2}(n) = (g_const+1) * alpha(2) * n(1);
% lambda_cell{3}(n) = g_const * alpha(3) * n(2);
% lambda_cell{4}(n) = alpha(4) * n(1);
% lambda_cell{5}(n) = alpha(5) * n(1) *(n(1)-1);
% lambda_cell{6}(n) = alpha(6) * n(1) *n(2);
% assumeAlso(alpha(1) == alpha(4));
% assumeAlso(alpha(2) == alpha(5));
% assumeAlso(alpha(3) == alpha(6));


Y = zeros(d, 2*K);
for j = 1:K
    Y(:, 2*j-1) = sources(:,j);
    Y(:, 2*j) = products(:,j);
end

YY = Y; % Copy the complex matrix

% Y = [0,0,0;1,0,0;2,0,0;1,1,0;0,1,1;0,0,1];  % Cycle
% Y = Y';

%% Find translated network without single complex merging.

Solution = {};
Index = {};

% [Solution,Index] = mergingcx(Y);  % Merging reactions
[Solution,Index] = CRN_translation(sources, products, 2);  % Merging reactions

% Sort out unique rows of Solution and Index
if numel(Solution) > 0
    [Solution,Index] = find_unique(Solution,Index);
end

% Solution{:}
% Index{:}

%% Find translated network considering single complex merging.

if numel(Solution) == 0
    for ii = 1:K
        ind1(ii) = {[ii]};
    end
    Y1 = YY;
    [net,ind] = singlecx(Y1,ind1);  % Merging a single cx and translate reactions
    % (e.g. when 2A->A+B & 0->A exist, 2A is merged into A
    % and reaction 2A->A+B is translated into A->B).
    Solution = [Solution;net];
    Index = [Index;ind];
else
    S = numel(Solution);
    for i = 1:S
        Y1 = cell2mat(Solution(i));
        ind1 = Index(i,:);
        [net,ind] = singlecx(Y1,ind1);   % Merging a single cx and translate reactions
        % (e.g. when 2A->A+B & 0->A exist, 2A is merged into A
        % and reaction 2A->A+B is translated into A->B).
        Solution = [Solution;net];
        Index = [Index;ind];
    end
end

% Sort out unique rows of Solution and Index
if numel(Solution) > 0
    [Solution,Index] = find_unique(Solution,Index);
end

%% Performing propensity factorization 
% Choose one of translated network
trans_net = 1; % index of translated network.
% 1 <= trans_net <= numel(Solution)

K_trans = sum(~cellfun(@isempty, Index(trans_net,:)));
% K_trans: the number of the reactions of the translated network
sol_Y = Solution{1};
sources_trans = sol_Y(:,1:2:(2*K_trans));
products_trans = sol_Y(:,2:2:(2*K_trans));

% construct propensities of the translated network.

ncell = sym2cell(n);
lambda_trans(n) = sym(zeros(K_trans,1));
lambda_trans_cell = sym2cell(formula(lambda_trans));

for k = 1:K_trans
    for j = cell2mat(Index(trans_net,k))
        lambda_trans_cell{k} = lambda_trans_cell{k} + lambda_cell{j};
    end
end

% set up the kinetic parameters kappa_k for the deterministic mass-action
% kinetic model on the translated network.

syms kappa [1 K_trans] positive
for k = 1:K_trans
    tmp = cell2mat(Index(trans_net,k));
    kappa(k) = alpha(tmp(end));
end

[a_list_tmp, b_list_tmp, elementary_basis] = CRN_find_elemtary_path(sources_trans); % the row vectors of H form a basis.
F_tmp = CRN_find_elementary_function(sources_trans, lambda_trans_cell, kappa);
for jj = 1:numel(F_tmp)
    F{jj} = simplify(F_tmp{jj});
end

% syms T0 positive
% assumeAlso(T0, 'integer')
start_point = ones(1,d);

elementary_coordinates = CRN_solve_sym_linear(elementary_basis, start_point);
coord_struct = struct2cell(elementary_coordinates);
theta(ncell{:}) = simplify(CRN_theta_construction(start_point, coord_struct, elementary_basis, F));
factorization_TF = CRN_check_factorization_condition(sources_trans, lambda_trans_cell, theta, kappa);


%% Compute CBE and derive a stationary distribuiton \pi(n).

[complexes_for_cbe, ~ , sources_idx] = unique(sources_trans, 'rows', 'stable');

num_C_trans = size(complexes_for_cbe,1);
products_idx = zeros(K_trans, 1);
for k = 1:K_trans
    for j = 1:num_C_trans
        if isequal(products_trans(k,:), complexes_for_cbe(j,:))
            products_idx(k) = j;
        end
    end
end

M_kappa = sym(zeros(num_C_trans, num_C_trans));

for k = 1:K_trans
   M_kappa(products_idx(k), sources_idx(k)) = kappa(k);
end

cbe = CRN_compute_cbe(complexes_for_cbe, M_kappa);

pi = simplify(prod(struct2cell(cbe)'.^n) / theta);


