function [d]=defi(Y)

N=numel(Y(1,:))/2; %%%%number of source complexes
M=numel(Y(:,1)); %%%number of species
Stoi=zeros(M,N);

for i=1:N
    Stoi(:,i)=Y(:,2*i)-Y(:,2*i-1);
end

C=unique(Y','rows')';
c=numel(C(1,:));

           
[S1,S2]=countlinkage(Y);
d=c-S2-rank(Stoi);

end

