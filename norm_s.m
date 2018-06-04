function S=norm_s(S,pop_s)
    m=max(pop_s);
    maxv=zeros(1,m);
    [n,~]=size(S);
    for i=1:n
        ind=pop_s(i);
        if maxv(ind)<S(i,1)
            maxv(ind)=S(i,1);
        end
    end
     for i=1:n
        ind=pop_s(i);
        S(i,1)=S(i,1)/maxv(ind);
     end
end