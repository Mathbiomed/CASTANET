function [X,Y]=find_unique(solution,index)
K=numel(solution);
smarker=[];
s=1;


%%%%%%%%Check which entry j of solution is equalt to entry k of solution
%%%%%%%% after an appropriate shuffle.


for k=1:K
    for j=k+1:K
        CS1=cell2mat(solution(k));  %%%% CS1 and CS2 are the complex matrices after all possible merging
        CS2=cell2mat(solution(j));     
        if numel(CS1(1,:))==numel(CS2(1,:))
            M=numel(CS1(1,:))/2;
            Shuffle3=perms(1:M);
            for i=1:numel(Shuffle3(:,1))          
                 for r=1:M                    
                    CS22(:,2*r-1)=CS2(:,2*Shuffle3(i,r)-1);  %%%clon of CS2 after shuffling
                    CS22(:,2*r)=CS2(:,2*Shuffle3(i,r));                       
                 end
                 if isequal(CS22,CS1)==1
                 smarker(s)=j;
                 s=s+1;
                 end
            end
        end      
    end
end
smarker=unique(smarker);



solution(smarker)=[];
index(smarker,:)=[];
X=solution;
Y=index;

end