function [d] = defi(sources, products)
Stoi = products - sources;
% the matrix whose column vectors are stoichiometric vectors of reactions.
C = unique([sources, products]','rows')';
c = size(C,2);
[~, S2] = countlinkage(sources, products);
d = c - S2 - rank(Stoi);
end