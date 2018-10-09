algorithm=@M2P_P;   
proname='DTLZ';
for q=1:6
    for M=[3,5,8,10]
        for t=[1,2]
            main('-algorithm',algorithm,'-problem',str2func([proname num2str(q)])...
                ,'-N',100,'-M',M,'-evaluation',100000,'-mode',2,'-run',t,'-M2P_P_parameter',10);
        end
    end
end
proname='WFG';
for q=1:9
    for M=[3,5,8,10]
        for t=[1,2]
            main('-algorithm',algorithm,'-problem',str2func([proname num2str(q)])...
                ,'-N',100,'-M',M,'-evaluation',100000,'-mode',2,'-run',t,'-M2P_P_parameter',10);
        end
    end
end 
obj=[3,5,8,10];
    proname='SDTLZ';
    for q=[1,3]
        obj_index = 1;
        for a=[10,5,3,2]
            M = obj(obj_index); 
            obj_index = obj_index + 1;
            for t = [1,2]
                main('-algorithm',@M2P_P,'-problem',str2func([proname num2str(q)])...
                       ,'-N',100,'-M',M,'-evaluation',100000,'-mode',2,'-run',t,'-M2P_P_parameter',{10,a}); 
            end
        end
    end      

