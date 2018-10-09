function pp=PP_value(S,wg_s,pop_s)
    global pro;
    S=norm_s(S,pop_s);
    [n,m]=size(S);
    up=zeros(n,m-1);
    for i=1:n
       ind=pop_s(i);
        for j=1:m-1
           d=S(i,j+1)-wg_s(ind,j+1);
           if d<0
               d=-d;
           end
           if d>up(ind,j)
                up(ind,j)=d;
           end
        end
    end
    pp=zeros(n,1);
    for i=1:n
        ind=pop_s(i);
        for j=1:m-1
            d=S(i,j+1)-wg_s(ind,j+1);
            if d<0
                d=-d;
            end
            pp(i)=pp(i)+2*sin(up(ind,j));
            pp(i)=pp(i)-sin(up(ind,j)-d);
            pp(i)=pp(i)-sin(d+up(ind,j));
        end
        pp(i)=(pp(i)+10e-6)*S(i,1);
        pp(i)=1/pp(i);
    end
end
