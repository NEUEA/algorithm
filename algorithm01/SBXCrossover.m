function child = SBXCrossover(parent1,parent2,pro)
%SBXCROSSOVER 此处显示有关此函数的摘要
%   此处显示详细说明
pd=pro.pd;
CrossProbability=1.0;
upperBound=pro.lb;
lowerBound=pro.ub;
distributionIndex=30;
EPS=1e-14;
child=zeros(2,pd);
if(rand(1)<=CrossProbability)
    for i=1:pd
        valueX1=parent1(i);
        valueX2=parent2(i);
        
        if rand(1)<=0.5 & valueX1-valueX2>EPS
            if valueX1<valueX2
                y1=valueX1;
                y2=valueX2;
            else
                y1=valueX2;
                y2=valueX1;
            end
            
            rnd=rand(1);
            beta=1.0+(2.0*(y1-lowerBound(i))/(y2-y1));
            alpha=2.0-(beta^(-(distributionIndex+1.0)));
            
            if(rnd<=(1.0/alpha))
                betaq=(rnd*alpha)^(1.0/(distributionIndex + 1.0));
            else
                betaq=(1.0 / (2.0 - rnd * alpha))^(1.0/(distributionIndex + 1.0));
            end
            c1=0.5*(y1+y2-betaq*(y2-y1));
            
            beta = 1.0 + (2.0 * (upperBound(i) - y2) / (y2 - y1));
            alpha=2.0-(beta^(-(distributionIndex + 1.0)));
            
            if(rnd<=(1.0/alpha))
                betaq=(rnd*alpha)^(1.0/(distributionIndex + 1.0));
            else
                betaq=(1.0 / (2.0 - rnd * alpha))^(1.0/(distributionIndex + 1.0));
            end
            c2 = 0.5 * (y1 + y2 + betaq * (y2 - y1));
            
            if(rand(1)<=0.5)
                child(1,i)=c2;
                child(2,i)=c1;
            else
                child(1,i)=c1;
                child(2,i)=c2;
            end
        else
            child(1,i)=valueX1;
            child(2,i)=valueX2;
        end
        
        for j=1:2
            if child(j,i) > upperBound(i)                     % 分量大于上限，设置为上限值
                child(j,i) = upperBound(i);
            elseif child(j,i) < lowerBound(i)                 % 分量小于下限，设置为下限值
                child(j,i) = lowerBound(i);
            end
        end
        
    end
end
end


