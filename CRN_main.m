% This is the main function for deriving stationary distribution. 
% The following four steps (Fig. 1) are sequentially performed by this code.
% 1) Network translation,
% 2) Propensity factorization,
% 3) CBE calculation
% 4) Deriving Stationary distribution.

%% Initialization of all input variables and parameters
clear; clc;

% Fig. 1 example
sources = [0 0; 1 0; 0 1; 2 0]';
products = [1 0; 0 1; 0 0; 1 1]';

% (number of reaction) * (number of species) matrix containing the source 
% and product complex vectors of reactions, respectively.

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

% New example 2  - cycle
% sources = [0 0 0; 1 1 0; 0 1 1]';
% products= [1 0 0; 0 2 0; 0 0 1]';



[d, K] = size(sources);
%propensities = [];
% source complexes of the reactions. (each complex could be repeated.)
% K and d could be determined by [K, d] = size(complexes) after determining
% complexes but for the simplicity user would enter them here, manually.

% d: the number of species
% K: the number of reactions.
syms n [d 1] integer 
% syms gg positive
assumeAlso(n >= 0) % nonnegativity of numbers
ncell = sym2cell(n);
lambda(n) = sym(zeros(K,1));
lambda_cell = sym2cell(formula(lambda));
% tt = 1;
syms alpha [K 1] positive

% Fig. 1 example
lambda_cell{1}(n) = alpha(1);
lambda_cell{2}(n) = alpha(2) * n(1);
lambda_cell{3}(n) = alpha(3) * n(2);
lambda_cell{4}(n) = alpha(4) * n(1) * (n(1) -1);

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

%% Performing Network translation

[Solution,Index] = CRN_translation(sources, products, 2);  % Merging reactions

disp(['The number of weakly reversible and deficiency zero translated networks is ', num2str(numel(Solution)/2),'.']);

%% Performing Propensity factorization 

% Choose one of translated network
PF_idx = 0;
for trans_net = 1:(numel(Solution)/2) % index of translated network.
    % 1 <= trans_net <= numel(Solution)
    disp(['Performing propensity factorization for a translated network ', num2str(trans_net), '.']);
    sources_trans = Solution{trans_net, 1};
    products_trans = Solution{trans_net, 2};
    stoi_trans = products_trans - sources_trans;

    K_trans = size(sources_trans, 2);
    % K_trans: the number of the reactions of the translated network

    % construct propensities of the translated network.

    % ncell = sym2cell(n);
    lambda_trans(n) = sym(zeros(K_trans,1));
    lambda_trans_cell = sym2cell(formula(lambda_trans));

    for k = 1:K_trans
        for j = Index{trans_net,k}
            lambda_trans_cell{k} = lambda_trans_cell{k} + lambda_cell{j};
        end
    end

    % set up the kinetic parameters kappa_k for the deterministic mass-action
    % kinetic model on the translated network.


    syms kappa [1 K_trans] positive
    for k = 1:K_trans
        tmp = cell2mat(Index(trans_net,k));
        kappa(k) = alpha(min(tmp));
    end

    [a_list_tmp, b_list_tmp, elementary_basis] = CRN_find_elementary_path(sources_trans); % the row vectors of H form a basis.
    F_tmp = CRN_find_elementary_function(sources_trans, lambda_trans_cell, ncell, kappa);
    for jj = 1:numel(F_tmp)
        F{jj} = simplify(F_tmp{jj});
    end

    % syms T0 positive
    % assumeAlso(T0, 'integer')
    sp_flag = 0; % flag for finding an appropriate start point n_0 for the theta construction step.
    sp_weight = 0;
    while sp_flag == 0
        sp_flag = 1;
        start_point = sp_weight * ones(d,1);
        for k = 1:K_trans
            sp_plus_stoi = num2cell(start_point + sources_trans(:,k));
            if ~isAlways(lambda_trans_cell{k}(sp_plus_stoi{:}) > 0)
                sp_flag = 0;
            end
        end
    end
    
    elementary_coordinates = CRN_solve_sym_linear(elementary_basis, start_point, n);

    conservation_law = null(stoi_trans', 'r')' * (n-start_point) == 0;
    assumeAlso(conservation_law)

    coord_cell = struct2cell(elementary_coordinates);
    theta(ncell{:}) = simplify(CRN_theta_construction(start_point, coord_cell, elementary_basis, F));
    factorization_TF = CRN_check_factorization_condition(sources_trans, lambda_trans_cell, theta, kappa);
    if prod(factorization_TF) == 1
        disp("All the factorization conditions hold!");
        PF_idx = 1;
        break;
    end
end

if PF_idx == 0
     disp("All the translated networks do not have the desired propensity factorization!");
end

%% Compute CBE and derive a stationary distribuiton pi(n).
[complexes_for_cbe, ~ , sources_idx] = unique(sources_trans', 'rows', 'stable');

complexes_for_cbe = complexes_for_cbe';

num_C_trans = size(complexes_for_cbe, 2);
products_idx = zeros(K_trans, 1);
for k = 1:K_trans
    for j = 1:num_C_trans
        if isequal(products_trans(:,k), complexes_for_cbe(:,j))
            products_idx(k) = j;
        end
    end
end

M_kappa = sym(zeros(num_C_trans, num_C_trans)); % 

for k = 1:K_trans
   M_kappa(products_idx(k), sources_idx(k)) = kappa(k);
end

cbe = CRN_compute_cbe(complexes_for_cbe, M_kappa);

pi = simplify(prod(struct2cell(cbe).^n) / theta);


disp("The source complexes of the translated network: ")
disp(Solution{trans_net,1})
disp("The product complexes of the translated network: ")
disp(Solution{trans_net,2})
for k = 1:K_trans
    if rem(k,10) == 1 & k ~= 11
        disp([num2str(k), 'st reaction propensity of the translated network: '])
    elseif rem(k,10) == 2 & k ~= 12
        disp([num2str(k), 'nd reaction propensity of the translated network: '])
    elseif rem(k,10) == 3 & k ~= 13
        disp([num2str(k), 'rd reaction propensity of the translated network: '])
    else
        disp([num2str(k), 'th reaction propensity of the translated network: '])
    end
    disp(lambda_trans_cell{k})
end
disp("Analytic formula of the stationary distribution: ")
disp(pi)

