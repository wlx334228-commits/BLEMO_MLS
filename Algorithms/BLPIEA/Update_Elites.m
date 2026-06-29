function  [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,W,T,beta,Maxgen,mode)
    
    Global = GLOBAL.GetObj();
    N = Global.N(3);
    switch mode
        case 1
            if ~isempty(Elite_Population)
                [Dis,inds] = min(pdist2(Elite_Population.upper_decs,History_data.upper_decs),[],2);
                Elite_Population(Dis==0) = [];
                Elite_Population = cat(2,Elite_Population,History_data(inds(Dis==0)).bilevel_best);
            else
                Elite_Population = History_data.bilevel_best;
            end
        case 2
            [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,W,T,beta,Maxgen,1);
            Elite_Population_temp = History_data.bilevel_best;
            if length(Elite_Population_temp) > Global.N(1)
                [~,Next] = Upper_Select(Elite_Population_temp,Global.N(1),W,'upper');
                [~,History_data] = lower_level_optimize(Elite_Population_temp(Next),History_data,T,beta,Maxgen);
            else
                [~,History_data] = lower_level_optimize(Elite_Population_temp,History_data,T,beta,Maxgen);
            end
            Elite_Population_temp = History_data.bilevel_best;
            Dis = min(pdist2(Elite_Population.upper_decs,Elite_Population_temp.upper_decs),[],2);
            Elite_Population(Dis==0) = [];
            Elite_Population = cat(2,Elite_Population,Elite_Population_temp);
    end
    
    if length(Elite_Population) > N
        Elite_Population = Upper_Select(Elite_Population,N,W,'bilevel');
    end
    
end
