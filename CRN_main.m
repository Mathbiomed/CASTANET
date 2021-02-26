% This is the main function for deriving stationary distribution. 
% The following four steps (Fig. 1) are sequentially performed by this code.
% 1) Network translation,
% 2) Propensity factorization,
% 3) CBE calculation
% 4) Deriving Stationary distribution.

%% Initialization of all input variables and parameters
clear; clc;

% Fig. 1 example
% sources = [0 0; 1 0; 0 1; 2 0]';
% products = [1 0; 0 1; 0 0; 1 1]';

% autophospho 3rd example
% sources = [1 0; 0 1; 1 1]';
% products = [0 1; 1 0; 0 2]';
% sources = [1 ; 0 ; 1 ]';
% products = [0 ; 1 ; 0 ]';

% example without suff condition
sources = [0 0; 1 0; 1 1; 1 0]';
products = [1 0; 1 1; 0 0; 2 0]';


% (number of reaction) * (number of species) matrix containing the source 
% and product complex vectors of reactions, respectively.

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
lambda_k = sym2cell(formula(lambda));
% tt = 1;
syms alpha [K 1] positive
% % Fig. 1 example
% lambda_k{1}(n) = alpha(1);
% lambda_k{2}(n) = alpha(2) * n(1);
% lambda_k{3}(n) = alpha(3) * n(2);
% lambda_k{4}(n) = alpha(4) * n(1) * (n(1) -1);
% 
% Autophospho 3rd example
% syms T0 integer
% assumeAlso(T0, 'positive');
% lambda_k{1}(n) = alpha(1) * n(1);
% lambda_k{2}(n) = alpha(2) * n(2);
% lambda_k{3}(n) = alpha(3) * n(1) * (T0 - n(1));

% example without suff condition
lambda_k{1}(n) = alpha(1);
lambda_k{2}(n) = alpha(2) * n(1);
lambda_k{3}(n) = alpha(3) * n(1) * n(2);
lambda_k{4}(n) = alpha(4) * n(1);


%% Performing Network translation

[Solution,Index] = CRN_translation(sources, products, 2);  % Merging reactions

disp(['The number of weakly reversible and deficiency zero translated networks is ', num2str(numel(Solution)/2),'.']);

%% Performing Propensity factorization 

% Choose one of translated network
PF_idx = 0;
for trans_net = 1:(numel(Solution)/2) % index of translated network.
    % 1 <= trans_net <= numel(Solution)
    disp(['Performing propensity factorization for the translated network ', num2str(trans_net), '.']);
    sources_trans = Solution{trans_net, 1};
    products_trans = Solution{trans_net, 2};
    stoi_trans = products_trans - sources_trans;

    K_trans = size(sources_trans, 2);
    % K_trans: the number of the reactions of the translated network

    % construct propensities of the translated network.

    % ncell = sym2cell(n);
    lambda_trans_tmp(n) = sym(zeros(K_trans,1));
    lambda_trans = sym2cell(formula(lambda_trans_tmp));

    for k = 1:K_trans
        for j = Index{trans_net,k}
            lambda_trans{k} = lambda_trans{k} + lambda_k{j};
        end
    end


    
    % set up the kinetic parameters kappa_k for the deterministic mass-action
    % kinetic model on the translated network.

    syms kappa [1 K_trans] positive
    
    for k = 1:K_trans
        try
            tmp = cell2mat(Index(trans_net,k));
            [coeff_tmp, term_tmp] = coeffs(lambda_trans{k}, n);
            coeff_tmp_list = coeff_tmp(ncell{:});
            coeff_of_highest_list = coeff_tmp_list(polynomialDegree(term_tmp) == max(polynomialDegree(term_tmp)));
            kappa(k) = coeff_of_highest_list(1);
            if isAlways(kappa(k) < 0)
                kappa(k) = -kappa(k);
            end
        catch
            kappa(k) = alpha(min(tmp));
        end
    end

    
    sp_flag = 0; % flag for finding an appropriate start point n_0 for the theta construction step.
    sp_weight = 0;
    while sp_flag == 0
        sp_flag = 1;
        start_point = sp_weight * ones(d,1);
        for k = 1:K_trans
            sp_plus_stoi = num2cell(start_point + sources_trans(:,k));
            if ~isAlways(lambda_trans{k}(sp_plus_stoi{:}) > 0)
                sp_flag = 0;
            end
        end
        if sp_weight > 100 % to avoid infinite while loop.
            disp('The weight of strat point is too large.');
            break
        end
        sp_weight = sp_weight + 1;
    end
        
    suff_cond_flag = 0; % flag for checking the sufficient condition for propensity factorization. 
    if prod(sum(sources_trans,1) <= 1)
        suff_cond_flag = 1;
        for k = 1:K_trans
            for ii = 1:d
                if isAlways(diff(lambda_trans{k}, n(ii)) == 0) & sources_trans(ii,k) == 1
                    suff_cond_flag = 0;
                elseif ~isAlways(diff(lambda_trans{k}, n(ii)) == 0) & sources_trans(ii,k) == 0
                    suff_cond_flag = 0;
                end
            end
        end     
    end
    
    if suff_cond_flag == 1   
        theta = sym(1);
        disp(['The sufficient condition for propensity factorization holds for the translated network ', num2str(trans_net), '.']);
        disp(['The factorization function theta(n) is constructed using the condition.']);
        disp(['The desired function theta is given by symprod(symprod(f_i(j),j,b_i+1,n_i), i, 1,', num2str(d), ') where']);
        syms j theta
        tmp = j*ones(d,1);
        tmp_cell = num2cell(tmp);
        for ii = 1:d
            disp(['b_',num2str(ii),' = ', num2str(start_point(ii))]);
        end
        disp('and');
        for ii = 1:d
            k = find(sources_trans(ii,:));
            disp(['f_',num2str(ii),'(j) = ']);
            disp(lambda_trans{k}(tmp_cell{:}));
        end
        PF_idx = 1;
        break;
    end
    
    clear j 
    
    
    % start_point = [5;5]; if stoi_trans == [];
    if prod(null(stoi_trans', 'r') > 0) == 1
        start_point = start_point + 5;
    end
    
        
    
    [a_list_tmp, b_list_tmp, elementary_basis] = CRN_find_elementary_path(sources_trans); % the row vectors of H form a basis.
    F_tmp = CRN_find_elementary_function(sources_trans, lambda_trans, ncell, kappa);
    for kk = 1:numel(F_tmp)
        F{kk} = simplify(F_tmp{kk});
    end

    elementary_coordinates = CRN_solve_sym_linear(elementary_basis, start_point, n);

    conservation_law = null(stoi_trans', 'r')' * (n-start_point) == 0;
    assumeAlso(conservation_law)
    if isa(elementary_coordinates, 'sym')
        coord_cell = sym2cell(elementary_coordinates);
    elseif isstruct(elementary_coordinates)
        coord_cell = struct2cell(elementary_coordinates);
    else
        disp('The type of elementary_coordinates is wrong.');
    end
    theta(ncell{:}) = simplify(CRN_theta_construction(start_point, coord_cell, elementary_basis, F));
    factorization_TF = CRN_check_factorization_condition(sources_trans, lambda_trans, theta, kappa, ncell);
    if prod(factorization_TF) == 1
        disp(['The factorization condition holds for the translated network ', num2str(trans_net), '!']);        
        PF_idx = 1;
        break;
    end
end

if PF_idx == 0
     disp("No translated network satisfying the factorization condition is identified.");
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
disp("The index of the translated network: ")

for k = 1:K_trans
    Index_trans{k} = Index{trans_net, k};
end

fprintf('      ');
for k = 1:K_trans
    fprintf('{');
    for jj = 1:length(Index_trans{k})
        tmplist = Index_trans{k};
        fprintf('%s',num2str(tmplist(jj)));
        if jj ~= length(Index_trans{k})
            fprintf(',');
        end
    end
    fprintf('} ');
end
fprintf('\n');

disp("ni is the number of the ith species.")
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
    disp(lambda_trans{k})
end
disp("Analytic formula of the stationary distribution: ")
disp(pi)

if suff_cond_flag == 1 
    disp('where theta is given above.');
end

