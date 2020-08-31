function [S1,S2] = countlinkage(Yout)
[N, M] = size(Yout);
% M: 2 * (the number of reactions);
% N: the number of species;
Y = [Yout;zeros(1,M)];

j = 1;
Y(end,1) = j;

% enumerate the complexes
for i = 2:M
    sw = 0;
    for k = 1:i-1
        if Y(1:N,i) == Y(1:N,k) & sw == 0
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