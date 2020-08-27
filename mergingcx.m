function [Solution,Index]=mergingcx(Y)



YY=Y;
N=numel(Y(1,:))/2; %%%%Number of reactions
shuffle=perms(1:N);
%M=numel(indexset(:,1)); %%%%Number of current groups

Solution={};
Index={};


s=1;


for m=1:numel(shuffle(:,1))
%%%%%%%%Start searching
Y=YY;
cxmarkers=[];
remarkers=[];
for ii=1:N
indexset(ii)={[ii]};
end


    for ii=1:N  %%%We fix a reaction to check whether merging is possible.  
         i=2*shuffle(m,ii)-1;     
          iindex=floor((i+1)/2);

          sw=0;         %%%searching switch


      for jj=ii+1:N  %%%Search a pair of complexes that have the same reaction vector as the fixed one.
          j=2*shuffle(m,jj)-1; 

      for kk=ii+1:N  
          k=2*shuffle(m,kk);      
          if isequal(Y(:,i+1)-Y(:,i),Y(:,k)-Y(:,j))==1   & (j~=i || k~=i+1)  & sw==0                 




              if k==j+1 %%% If the pair forms an existing reaction                  
                  sw=1;             
                  jindex=floor((j+1)/2);             
                  cxmarkers=[cxmarkers,i,i+1];  %%% mark which complexes are merged
                  remarkers=[remarkers,iindex];  %%% mark which reactions are merged
              %%%%%%%Update the indexset by adding additional merging index        

                  v1=cell2mat(indexset(iindex));
                  v2=cell2mat(indexset(jindex));             
                  indexset(jindex)={union(v1,v2)};             
              end



              %%%%%%%%%%Check the current translation gives us def0 and wr
                 if sw==1 
                   indexset_tem=indexset;
                   Y_tem=Y;
                   indexset_tem(remarkers)=[];
                   Y_tem(:,cxmarkers)=[];
                   [S1,S2]=countlinkage(Y_tem);
                if defi(Y_tem)==0 && S1==S2    

                      Solution=[Solution;Y_tem];
                      for jjj=1:numel(indexset_tem)
                      Index(s,jjj)={cell2mat(indexset_tem(jjj))};
                      end
                      s=s+1;
                  end
               end



          end     

      end
      end


    if sw==0  %%%%if no reaction merging has been occured, we check complex merging.

        for j=i+2:2*N  %%%Search a pair of complexes that have the same reaction vector as the fixed one.
          for k=i+2:2*N  
    %
               sw=0;
               if Y(:,i+1)-Y(:,i)==Y(:,k)-Y(:,j)  & sw==0  & isequal(Y(:,i),Y(:,j))==0              
                 I=cell2mat(indexset(iindex))  ; 

                  for ii=1:numel(I)
                  Y(:,ii+1)=Y(:,k);         
                  Y(:,ii)=Y(:,j);     
                  end
                  sw=1; %%%% sw=1 stops unnecessary additional merging for reaction i
               end




                 %%%%%%%%%%Check the current translation gives us def0 and wr
                  if sw==1 
                     %[i,j,k]
                    indexset_tem=indexset;
                    Y_tem=Y;
                    indexset_tem(remarkers)=[];
                    Y_tem(:,cxmarkers)=[];
                    [S1,S2]=countlinkage(Y_tem);
                   if defi(Y_tem)==0 && S1==S2    

                      Solution=[Solution;Y_tem];
                      for jjj=1:numel(indexset_tem)
                      Index(s,jjj)={cell2mat(indexset_tem(jjj))};
                      end
                      s=s+1;
                   end
                  end


          end
          end
    end





    end



end


end



