clear
% Y=[2,0,0;1,1,0;0,1,0;1,0,0;0,1,0;0,0,1;0,0,1;0,1,0];  %%% Figure 2d
% Y=Y';
%Y=[2,0;1,1;1,1;1,0;1,0;2,0]         

Y=[2,0;1,1;0,0;1,0;1,0;0,1;0,1;0,0];    %%%% Figure 1
Y=Y';
% 
% Y=[0,0,0;1,0,0;2,0,0;1,1,0;0,1,1;0,0,1];  %%%Cycle
% Y=Y';



YY=Y; %%%%%%%%%Copy the complex matrix
N=numel(Y(1,:))/2; %%%% Number of source complexes
Solution={};
Index={};





    
[Solution,Index]=mergingcx(Y);  %%%%%%Merging reactions




%%%%%%%% Sort out unique rows of Solution and Index
if numel(Solution)>0
[Solution,Index]=find_unique(Solution,Index);  
end

Solution

%%%%%%%%%%% Single complex merging
if numel(Solution)==0          
    for ii=1:N
    ind1(ii)={[ii]};
    end  
    Y1=YY;    
    [net,ind]=singlecx(Y1,ind1);  %%%%Merging a single cx and translate reactions 
                                     %%%% (e.g. when 2A->A+B & 0->A exist, 2A is merged into A
                                     %%%% and reaction 2A->A+B is translated into A->B). 
    
    Solution=[Solution;net];
    Index=[Index;ind];   
else
        S=numel(Solution);
        for i=1:S
            Y1=cell2mat(Solution(i));
            ind1=Index(i,:);
            [net,ind]=singlecx(Y1,ind1);   %%%%Merging a single cx and translate reactions 
                                             %%%% (e.g. when 2A->A+B & 0->A exist, 2A is merged into A
                                             %%%% and reaction 2A->A+B is translated into A->B). 
       
        Solution=[Solution;net];
        Index=[Index;ind];      
        end
end





%%%%%%%% Sort out unique rows of Solution and Index
if numel(Solution)>0
[Solution,Index]=find_unique(Solution,Index);
end












    