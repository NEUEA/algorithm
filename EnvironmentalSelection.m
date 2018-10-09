function [Population,MMF,p] = EnvironmentalSelection(Population,N,p,M)
% The environmental selection of MD-MOEA

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    MMF=[M2F_p(Population,p),(1:size(Population,2))'];
    MMF=sortrows(MMF,1);
    limit=10^-2;    
    if(size(find(MMF(:,1)<(-(10^-5))),1)<N)
        p=p-0.1
    else 
        p=p+0.1
    end
    if(p<0)
        p=0;
    end
    s=1; i=1; P=[];
    non_Pop=MMF(find(MMF<0),2);
    select=zeros(1,N);
    %limit=min(abs(MMF(:,1)))*10^-3;
    while((s<=N)&&(i<=size(Population,2)))
        if is_similarity(Population,i,select,limit)
            select(s)=MMF(i,2);
            s=s+1;
            if(i<=size(non_Pop,1))
                non_Pop(i)=0;
            end
            
        end
        i=i+1;
    end
    MMF(select(1:s-1),:)=[];
    if s<=N
         i=1;
         while s<=N
             select(s)=MMF(i,2);
             s=s+1;
             i=i+1; 
         end
    else
        not_select_Pop=non_Pop(find(non_Pop>0));
        if(size(not_select_Pop,1)>0)
            for i=1:size(not_select_Pop,1)
                if is_similarity(Population,not_select_Pop(i),select,limit)
                    [dItoX,X]=nearest(Population(not_select_Pop(i)).obj,Population(select).objs,0);
                    R=randperm(size((select),2),1);
                    while R==X
                        R=randperm(size((select),2),1);
                    end
                    [dRtoY,Y]=nearest(Population(select(R)).obj,Population(select).objs,R);
                    if dItoX>dRtoY
                        select(R)=not_select_Pop(i);
                    else
                        [dXtoZ,Z]=nearest(Population(select(X)).obj,Population(select).objs,X);
                        [dItoW,W]=nearest(Population(not_select_Pop(i)).obj,Population(select).objs,X);
                        if dItoW>dXtoZ
                            select(X)=not_select_Pop(i);
                        end
                    end
                end
            end
        end
        
    end
    Population=Population(select);
    MMF=M2F_p(Population,p)';

end