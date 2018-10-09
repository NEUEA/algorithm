function f = M2F_p( Population,p )
%M2F_P 此处显示有关此函数的摘要
%   此处显示详细说明

max_x=zeros(length(Population),1);
objs=Normalization(Population.objs);
%objs=Population.objs;
for i=1:size(Population,2)
    x1=repmat(objs(i,:),length(Population),1);
    min_x=min(x1-objs,[],2).*((abs(min(x1-objs,[],2))./sum(abs(x1-objs),2)).^p);
    max_x(i)=max(min_x);
end
   f=max_x;
end


