function [solution,index]=singlecx(Y,ind)
[M, tmp] = size(Y);
N = tmp/2;
% M: the number of species.
% N: the number of reactions

solution={};
index={};

Shuffle2=perms(1:2*N);
% Now we check a single complex merging

YY=Y;
for k=1:numel(Shuffle2(:,1))

   Y=YY; % Set Y to be the original Y. 
   Merging_check=zeros(2*N,1); % i & i+1 th entry of this array are 1 once either i th cx or i+1 th cx was merged to another cx. 
    for ii=1:2*N
        i=Shuffle2(k,ii);  % fix a complex
        for jj=1:2*N % Search a complex that allows a reaction shifting by mering a single complex 
            j=Shuffle2(k,jj); % searching cx through the shuffled enumeration of cx'                 
            if j~=i
                eta=Y(:,j)-Y(:,i); 
                 if mod(i,2)==1  % If the i the complex is a source      
                       if Y(:,i+1)+eta>=0 && Merging_check(i+1,1)==0                 

                           % Translate cx i and cx i+1

                              Y(:,i+1)=Y(:,i+1)+eta;         
                              Y(:,i)=Y(:,j);     
                              Merging_check(i)=1;
                              Merging_check(i+1)=1;
                        end
                 else          % If the i the complex is a product
                        if Y(:,i-1)+eta>=0 && Merging_check(i-1,1)==0             
                           % Translate cx i and cx i+1
                              Y(:,i-1)=Y(:,i-1)+eta;         
                              Y(:,i)=Y(:,j);     
                              Merging_check(i)=1;
                              Merging_check(i-1)=1;
                        end            
                 end
              
                             % Check the current translation gives us def0 and wr
                      
                  [S1,S2]=countlinkage(Y);
                   
                   if defi(Y)==0 && S1==S2    
                      x=check_unique(Y,solution)     ;             
                      if x==1
                          123213
                      solution=[solution;{Y}]; 
                      index=[index;ind];                     
                      end
                   end               
            end
        end
    end
  
end

11111
Y=YY; % Refresh Y to the original Y. 
% Now we cancel out redundant species for each reaction.
for i=1:2:2*N

    for m=1:M
       t=Y(m,i)-Y(m,i+1);
        if t>=0
            Y(m,i)=t;
            Y(m,i+1)=0;
        else
            Y(m,i)=0;
            Y(m,i+1)=-t;
        end
    end

    % Check the current translation gives us def0 and wr          

                [S1,S2]=countlinkage(Y);
               if defi(Y)==0 && S1==S2    
                  x=check_unique(Y,solution); 
                  if x==1
                  solution=[solution;Y];  
                  index=[index;ind];                     
                  end
               end
end

end



