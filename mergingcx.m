function [Solution,Index] = mergingcx(Y)

K = size(Y,2)/2; % the number of reactions
d = size(Y,1); % the number of species (dimension)
sources_origin = Y(:, 1:2:(2*K));
products_origin = Y(:, 2:2:(2*K));

Y_origin = Y;
shuffle=perms(1:K);

Solution = {};
Index = {};
s = 1;

for m = 1:numel(shuffle(:,1))
% for m = 1:15
    % Start searching
    sources = sources_origin;
    products = products_origin;
    remarkers = [];
    for ii = 1:K
        indexset(ii) = {[ii]};
    end
    
    for ii = 1:K  % We fix a reaction to check whether merging is possible.
        i = shuffle(m,ii);       
        sw = 0;% searching switch
        for jj = (ii+1):K  % Search a pair of complexes that have the same reaction vector as the fixed one.
            j = shuffle(m,jj);
            if i > size(sources, 2) | j > size(sources, 2)
                continue
            end
            if isequal(products(:,i)-sources(:,i), products(:,j)-sources(:,j))% & sw == 0
                % merging ith reaction to jth reaction
                remarkers = i;  % mark which reactions are merged
                % Update the indexset by adding additional merging index
                
                v1 = cell2mat(indexset(i));
                v2 = cell2mat(indexset(j));
                indexset(j)={[v1, v2]};

                sources(:, i) = sources(:, j);
                products(:, i) = products(:, j);
                % Check the current translation gives us def0 and wr
                indexset(remarkers) = [];
                sources(:, remarkers) = [];
                products(:, remarkers) = [];
                K_trans = size(sources, 2);
                Y_tmp = zeros(d, 2*K_trans);
                for idx = 1:K_trans
                    Y_tmp(:, 2*idx-1) = sources(:,idx);
                    Y_tmp(:, 2*idx) = products(:,idx);
                end
                [S1,S2] = countlinkage(sources, products);
                if defi(sources, products) == 0 & S1 == S2
                    Solution = [Solution;Y_tmp];
                    for jjj = 1:numel(indexset)
                        Index(s,jjj) = {cell2mat(indexset(jjj))};
                    end
                    s = s + 1;
                end
                
            end
        end
        if sw == 0  % if no reaction merging has been occured, we check complex merging.
            for j = 1:(2*size(sources,2))  % Search a pair of complexes that have the same reaction vector as the fixed one.
                for k = 1:(2*size(sources,2))
                    sw = 0;
                    complexes_tmp = [sources, products];
                    if i > size(sources,2) % | i > numel(indexset)
                        continue
                    end
                    if isequal(products(:,i) - sources(:,i), complexes_tmp(:,k) - complexes_tmp(:,j)) & ~isequal(products(:,i), complexes_tmp(:,k)) & k-j ~= K
                        products(:, i) = complexes_tmp(:,k);
                        sources(:, i) = complexes_tmp(:,j);
                        % Check the current translation gives us def0 and wr

                        K_trans = size(sources,2);
                        Y_tmp = zeros(d, 2*K_trans);
                        for idx = 1:K_trans
                            Y_tmp(:, 2*idx-1) = sources(:,idx);
                            Y_tmp(:, 2*idx) = products(:,idx);
                        end
                        [S1,S2] = countlinkage(sources, products);
                        if defi(sources, products) == 0 & S1 == S2
                            Solution = [Solution;Y_tmp];
                            for jjj = 1:numel(indexset)
                                Index(s,jjj) = {cell2mat(indexset(jjj))};
                            end
                            s = s + 1;
                        end
                    end
                end
                sw = 1; %  sw = 1 stops unnecessary additional merging for reaction i
            end
        end
  
    end
    
end

end
