function [x]=check_unique(Y,solution)
K=numel(solution);



%%%%%%%%Check which entry j of solution is equalt to entry k of solution
%%%%%%%% after an appropriate shuffle.
sw=0;
CS1=Y;  %%%% CS1 and CS2 are the complex matrices after all possible merging
    for j=1:K
        
        CS2=cell2mat(solution(j));     
        if numel(CS1(1,:))==numel(CS2(1,:)) && sw==0
            M=numel(CS1(1,:))/2;
            Shuffle3=perms(1:M);
            for i=1:numel(Shuffle3(:,1))          
                 for r=1:M                    
                    CS22(:,2*r-1)=CS2(:,2*Shuffle3(i,r)-1);  %%%clon of CS2 after shuffling
                    CS22(:,2*r)=CS2(:,2*Shuffle3(i,r));                       
                 end
                 if isequal(CS22,CS1)==1                
                 sw=1;
                 end
            end
        end      
    end
    
    if sw==1   %%%%Solution is duplicated
        x=0;
    else 
       x=1;
    end
    
end






