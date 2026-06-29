function Offspring = Generate(Population_B,History_data,redundant_dec,W,T)
    Global = GLOBAL.GetObj();
    
    redundant_dec_u = redundant_dec(redundant_dec<=Global.D(1));
    redundant_dec_l = redundant_dec(redundant_dec>Global.D(1))-Global.D(1);
    
    N = length(Population_B);
    
%     History_Population = NSGAII_Update(History_data.Best_Archives,N,'upper');
    
    [~,History_Population] = Update_Data(History_data,N);
    
    Distance = min(pdist2(Population_B.upper_decs,History_Population.upper_decs),[],2);
    
    MatingPop = [Population_B(Distance>0),History_Population];
    
    Target = Upper_Mating(MatingPop,N,W,T,2);
    
    Pool = PoolMating(Target,History_Population,T,2);
    Target_decs = Target.upper_decs;
    upper_decs = DE_current_1_bin(Target_decs,Pool(:,1).upper_decs,Pool(:,2).upper_decs,'upper',{0.5,0.5,1,20});
    upper_decs(:,redundant_dec_u) = Target_decs(:,redundant_dec_u);
    
%     [~,LinkID] = min(pdist2(s,History_Population.upper_decs),[],2);
%     Target = Global.Evaluate(INDIVIDUAL({s,History_Population(LinkID).lower_decs}),'upper');
%     
%     Pool = PoolMating(Target,History_Population,T,2);
    Target_decs = Target.lower_decs;
    lower_decs = DE_current_1_bin(Target.lower_decs,Pool(:,1).lower_decs,Pool(:,2).lower_decs,'lower',{0.5,0.2,0,20});
    lower_decs(:,redundant_dec_l) = Target_decs(:,redundant_dec_l);
    
    Offspring = INDIVIDUAL({upper_decs,lower_decs});
    Offspring = Global.Evaluate(Offspring,'upper');
    
end

function Pool = PoolMating(Target,MatingPop,T,num)
    N = length(Target);
    B = pdist2(Target.upper_objs,MatingPop.upper_objs);
    B(B==0) = inf;
    [~,B] = sort(B,2);
    B = B(:,1:min(T,length(MatingPop)));
    Pool = [];
    for i=1:N
        if rand < 0.5
            Pool = cat(1,Pool,B(i,randperm(T,num)));
        else
            Pool = cat(1,Pool,randperm(length(MatingPop),num));
        end
    end
    Pool=MatingPop(Pool);
end
