function [d,x]= nearest( i,P,a )
%NEAREST �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    if a>0
    P(a,:)=[];
    end
    P=bsxfun(@minus, P, i);
    d=sum(P.*P,2);
    [d,x]=min(d);
    if (x>=a&&a~=0)
        x=x+1;
    end
end

