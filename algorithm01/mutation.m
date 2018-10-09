function f=mutation(child,p,pro,mum)
low=pro.ub;
up=pro.lb;
if up==low
    f=child;
    return
end
for i=1:pro.pd
    if rand(1)<p
        y=child(i);
        delta1=(y-low(i))/(up(i)-low(i));
        delta2=(up(i)-y)/(up(i)-low(i));
        rnd=rand(1);
        mutpow=1/(mum+1);
        if rnd<=0.5
            xy=1-delta1;
            val=2*rnd+(1-2*rnd)*(xy^(mum+1));
            deltaq=val^mutpow-1;
        else
            xy=1-delta2;
            val=2*(1-rnd)+2*(rnd-0.5)*(xy^(mum+1));
            deltaq=1-val^mutpow;
        end
        y=y+deltaq*(up(i)-low(i));
        if y>up(i)
            y=up(i);
        end
        if y<low(i)
            y=low(i);
        end
        child(i)=y;
    end
    
end
f=child;
end
        
    
