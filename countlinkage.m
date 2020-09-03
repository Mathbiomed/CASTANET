function [S1,S2] = countlinkage(sources, products)
[d, K] = size(sources); 
%K: the number of reactions
%d: the number of species (dimension)
Yout = zeros(d, 2*K);
for idx = 1:K
    Yout(:, 2*idx-1) = sources(:,idx);
    Yout(:, 2*idx) = products(:,idx);
end
M = 2*K;

Y = [Yout;zeros(1,M)];

j = 1;
Y(end,1) = j;

% enumerate the complexes
for i = 2:M
    sw = 0;
    for k = 1:i-1
        if Y(1:d,i) == Y(1:d,k) & sw == 0
            Y(end,i) = Y(end,k);
            sw = 1;
        end
    end
    if sw == 0
        j = j+1;
        Y(end,i) = j;
    end
end

C = zeros(j,j);
for i = 1:2:M
    n = Y(end,i);
    m = Y(end,i+1);
    C(n,m) = 1;
end

A = sparse(C);
[S1, ~] = graphconncomp(A);               % number of strongly connect component
[S2, ~] = graphconncomp(A,'weak','true'); % number of linkage classes

end