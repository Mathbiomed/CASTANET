function [d] = defi(Y)

[M, tmp] = size(Y);
N = tmp/2;
% M: the number of species.
% N: the number of reactions == the number of source complexes

Stoi = zeros(M, N);

for i = 1:N
    Stoi(:,i) = Y(:,2*i) - Y(:,2*i-1);
end

C = unique(Y','rows')';
c = numel(C(1,:));

           
[S1,S2] = countlinkage(Y);
d = c - S2 - rank(Stoi);

end