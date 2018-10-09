function APIgEA(proname,od,time)
global pro;
global B_Angle;
global count;
count=1;
oldcount=count;
pro.name=proname;
pro.od=od;
maxCycle=1500;
digits(50)
switch pro.od
    case 3            % DTLZ1
        pro.H=12;
    case 5
        pro.H=5;
    case 8
        pro.H=3;
    case 10
        pro.H=3;
end
switch pro.od
    case 3
        pro.popsize=100;
        B_Angle=[0,0.451026811796262,0.643501108793284,0.795398830184144,0.927295218001612,1.04719755119660,1.15927948072741,1.26610367277950,1.36943840600457,1.47062890563334,1.57079632679490;
            0,0.157079632679490,0.314159265358979,0.471238898038469,0.628318530717959,0.785398163397448,0.942477796076938,1.09955742875643,1.25663706143592,1.41371669411541,1.57079632679490];
    case 5
        pro.popsize=256;
         B_Angle=[0,0.978973572693096,1.21610987017399,1.40173777042144,1.57079632679490;
            0,0.883170143467152,1.15494073013670,1.37184244880235,1.57079632679490;
            0,0.722734247813416,1.04719755119660,1.31811607165282,1.57079632679490;
            0,0.392699081698724,0.785398163397448,1.17809724509617,1.57079632679490];
    case 8
        %         pro.popsize=128;
        %         B_Angle=[0,1.30821604508933,1.57079632679490;
        %             0,1.28582735936785,1.57079632679490;
        %             0,1.25657868559176,1.57079632679490;
        %             0,1.21610987012030,1.57079632679490;
        %             0,1.15494072999723,1.57079632679490;
        %             0,1.04719755119660,1.57079632679490;
        %             0,0.785398163397448,1.57079632679490];
        pro.popsize=243;
        B_Angle=[0,1.12430643928177,1.36908524593775,1.57079632679490;
            0,1.06820050306718,1.34275118063434,1.57079632679490;
            0,0.984484356570793,1.30266283736669,1.57079632679490;
            0,0.841068670567930,1.23095941734077,1.57079632679490;
            0,0.523598775598299,1.04719755119660,1.57079632679490];
    case 10
          pro.popsize=256;
%         B_Angle=[0,1.12430643928177,1.36908524593775,1.57079632679490;
%             0,1.06820050306718,1.34275118063434,1.57079632679490;
%             0,0.984484356570793,1.30266283736669,1.57079632679490;
%             0,0.841068670567930,1.23095941734077,1.57079632679490;
%             0,0.523598775598299,1.04719755119660,1.57079632679490];
%        B_Angle=[0,1.34070403532739,1.57079632679490;
%             0,1.32605731853570,1.57079632679490;
%             0,1.30821604508933,1.57079632679490;
%             0,1.28582735936785,1.57079632679490;
%             0,1.25657868559176,1.57079632679490;
%             0,1.21610987012030,1.57079632679490;
%             0,1.15494072999723,1.57079632679490;
%             0,1.04719755119660,1.57079632679490;
%             0,0.785398163397448,1.57079632679490];
        B_Angle=[0,1.34070403532739,1.57079632679490;
            0,1.32605731853570,1.57079632679490;
            0,1.57079632679490,1.57079632679490;
            0,1.28582735936785,1.57079632679490;
            0,1.25657868559176,1.57079632679490;
            0,1.21610987012030,1.57079632679490;
            0,1.57079632679490,1.57079632679490;
            0,1.04719755119660,1.57079632679490;
            0,0.785398163397448,1.57079632679490];
end
% EF=load(['PF\' pro.name '-' num2str(pro.od) '-testPF.pf']);
EF=load(['PF\' pro.name '(' num2str(pro.od) ')-PF.txt']);
switch pro.name
    case 'DTLZ1'
        pro.pd=pro.od+5;
        pro.ub=zeros(1,pro.pd);
        pro.lb=ones(1,pro.pd);
        ref=ones(1,pro.od);
    case {'DTLZ2','DTLZ3','DTLZ4'}
        pro.pd=pro.od+10;
        pro.ub=zeros(1,pro.pd);
        pro.lb=ones(1,pro.pd);
        ref=ones(1,pro.od);
    otherwise %WFG
        pro.pd=pro.od*2-1+20;
        pro.ub=zeros(1,pro.pd);
        pro.lb=2:2:pro.pd*2;
        ref=(2:2:pro.od*2)+1;
        pro.k=(pro.od-1)*2;
end
mu=30;
mum=20;
for t=1:time
    wg=init_weights(pro);
    [pro.popsize,~]=size(wg);
    pop=initialize_variables(pro);
   
    iter=1;
    while iter<=maxCycle
        tic;
       
        parent_chromosome=mate_selection(pop);
        offspring=genetic_operator(parent_chromosome,pro,mu,mum);
        [n,m]=size(offspring);
        val=zeros(n,pro.od);
        for i = 1:n
            val(i,1:pro.od)=evaluate(offspring(i,:),pro);
        end
        offspring=[offspring,val];
        next_pop=[pop;offspring];
        
        S=normalization(next_pop);
        for i=1:pro.popsize
            wg_s(i,:)=cartesian_polar(wg(i,:),[]);
        end
        pop_s=Partition_Association(next_pop,wg);%每个基因属于哪个权重
        
        pop=zeros(pro.popsize,pro.pd+pro.od);
        %     ARC=ARC_selection(S,pro.popsize/2);
        %     pop(1:pro.popsize/2,:)=next_pop(ARC,:);
        %     S(ARC,:)=[];
        PP=PP_selection(S,wg_s,pop_s);
        pop(1:pro.popsize,:)=next_pop(PP,:);
        run_time(oldcount:count)=toc;
        i_dg(oldcount:count)=IGD(EF,pop(:,pro.pd+1:pro.pd+pro.od));
        hv(oldcount:count)=sum(hypeIndicatorSampled(pop(:,pro.pd+1:pro.pd+pro.od),ref,pro.popsize,10000))/prod(ref);
        oldcount=count;
        disp(sprintf('iteration %u finished, time used: %u HV %u IGD %u', iter, run_time(count),hv(count),i_dg(count)));
            plot(1:pro.od,pop(:,end-pro.od+1:end));
            drawnow
        iter=iter+1;
        
    end
    save(['测试结果\' pro.name '_' num2str(pro.od) '_' num2str(t) '.mat']);
     plot(1:pro.od,pop(:,end-pro.od+1:end));
            xlabel('object-number');
            ylabel('object-value');
            title(proname);
            print(gcf,'-djpeg',['img\' pro.name '_' num2str(pro.od) '_' num2str(t) '.jpg']);

end
clear global;
end
