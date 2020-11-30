function [Solution,Index] = CRN_translation(sources, products, max_order)
if nargin == 2
    max_order = 2;
end

[d, K] = size(sources);
% d: the number of species (dimension)
% K: the number of reactions

% if all the reactions only consume (or produce) certain species
% than the network cannot be weakly reversible even with network
% translation.
stoi_vec = products - sources;
is_consume = (stoi_vec <= 0);
is_produce = (stoi_vec >= 0);
is_preserve = (stoi_vec == 0);
stoi_dim = rank(stoi_vec);
if sum(prod(is_consume,2) + prod(is_produce,2) -2*prod(is_preserve,2)) > 0
   Solution = {};
   Index = {};
   return
end

sources_reduced = -ones(d, K);
products_reduced = -ones(d, K);

for k = 1:K
    sources_reduced(:, k) = sources(:,k) - min(sources(:,k), products(:,k));
    products_reduced(:, k) = products(:,k) - min(sources(:,k), products(:,k));
end

% For each reaction, sources_simplified and products_simplified have
% disjoint supports because all common species are subtracted.

tmp = nchoosek(1:(d+max_order),d);
max_copy_number = size(tmp,1);
additional_complex_tmp = (diff([zeros(max_copy_number,1), tmp], 1, 2) - 1)';
additional_complex = -ones(d, nchoosek(d+max_order-1,d-1), max_order);
order_number = zeros(max_order, 1);
for j = 2:size(additional_complex_tmp,2)
    order_tmp = sum(additional_complex_tmp(:,j));
    order_number(order_tmp) = order_number(order_tmp) + 1;
    additional_complex(:,order_number(order_tmp),order_tmp) = additional_complex_tmp(:,j);
end

% 
sources_copy = -ones(d, max_copy_number, K);
products_copy = -ones(d, max_copy_number, K);
order_number_cumul = cumsum(order_number);
% reaction_orders = sum(sources_reduced); % consider only sources for the reaction orders
reaction_orders = max(sum(sources_reduced), sum(products_reduced)); % consider both sources and products for the reaction orders. It is valid because the resulting network is weakly reversible.

for k = 1:K
    sources_copy(:,1,k) = sources_reduced(:,k);
    products_copy(:,1,k) = products_reduced(:,k);
    for m = 1:(max_order - reaction_orders(k))
        if m == 1
            sources_copy(:,2:(order_number_cumul(m)+1),k) = sources_reduced(:, k) + additional_complex(:,1:order_number(m),m);
            products_copy(:,2:(order_number_cumul(m)+1),k) = products_reduced(:, k) + additional_complex(:,1:order_number(m),m);
        else
            sources_copy(:,(order_number_cumul(m-1)+2):(order_number_cumul(m)+1),k) = sources_reduced(:, k)  + additional_complex(:,1:order_number(m),m);
            products_copy(:,(order_number_cumul(m-1)+2):(order_number_cumul(m)+1),k) = products_reduced(:, k)  + additional_complex(:,1:order_number(m),m);
        end
    end
end
number_of_copies = zeros(1,K);
for k = 1:K
    if reaction_orders(k) == max_order
        number_of_copies(k) = 1;
    else
        number_of_copies(k) = order_number_cumul(max_order - reaction_orders(k)) + 1;
    end
end
number_of_comb = prod(number_of_copies)

% reject if there is no possibility to attain the desired properties:
% weak reversibility and zero deficiency.




%
%

sources_tmp = -ones(d, K);
products_tmp = -ones(d, K);

sol_idx = 1;

Solution = {};
Index = {};

current_copy = ones(K, 1);
current_copy(1) = 0;
% create a matrix containing a vector between ones(K,1) and number_of_copies.
for ii = 1:number_of_comb
    if rem(ii,10000) == 0
        disp('ii');
        disp(ii);
    end
    current_copy(1) = current_copy(1) + 1;
    for k = 1:K
        if current_copy(k) > number_of_copies(k)
            current_copy(k) = 1;
            current_copy(k+1) = current_copy(k+1) + 1;
        end
    end 
    
    for k = 1:K
        sources_tmp(:, k) = sources_copy(:, current_copy(k), k);
        products_tmp(:, k) = products_copy(:, current_copy(k), k);
    end
    
    num_complexes = size(unique([sources_tmp, products_tmp]','rows')',2);
    
    % if any of the below conditions holds then the deficiency cannot be 0.
    if num_complexes < stoi_dim + 1
        continue
    elseif num_complexes > 2*stoi_dim
        continue
    end
    
    % the below condition indicate that there is unused species.
    if sum(prod([sources,products] == 0, 2)) > 0
        continue
    end 
    
    [S1,S2] = countlinkage(sources_tmp, products_tmp);
    % S1: the number of strongly connected components
    % S2: the number of linkage classes
    
    deficiency = num_complexes - S2 - stoi_dim;
    if deficiency == 0 & S1 == S2
        sol_idx = sol_idx + 1;
        complexes_tmp = [sources_tmp; products_tmp];
        [complexes_unique, ~, ic] = unique(complexes_tmp', 'rows');
        complexes_unique = complexes_unique';
        sources_unique = complexes_unique(1:d,:);
        products_unique = complexes_unique((d+1):(2*d),:);
        K_trans = size(sources_unique, 2);
        Solution = [Solution; {sources_unique, products_unique}];
        indexset = cell(1,K_trans);
        for k = 1:K
            indexset{ic(k)} = [indexset{ic(k)}, k];
        end
        Index = [Index; indexset];
    end
%     if rem(ii,floor(number_of_comb/20)) == 0
%          disp(['Searching translated netowrks ... ', num2str(100*ii/number_of_comb, 3), '% completed.']);
%     end
end

end
