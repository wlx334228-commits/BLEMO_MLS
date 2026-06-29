function  [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,redundant_dec_lobj,W,T,alpha,belta,Maxgen,mode)
    
    Global = GLOBAL.GetObj();
    N = Global.N(3);
    switch mode
        case 1
            if ~isempty(Elite_Population)
                [Dis,inds] = min(pdist2(Elite_Population.upper_decs,History_data.upper_decs),[],2);
                Elite_Population(Dis==0) = [];
                Elite_Population = [Elite_Population,History_data(unique(inds(Dis==0))).bilevel_best];
            else
                Elite_Population = History_data.bilevel_best;
            end
            
        case 2
            
            [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,redundant_dec_lobj,W,T,alpha,belta,Maxgen,1);
            
            [~,History_Population] = Update_Data(History_data,Global.N(1));
            Dis_temp = min(pdist2(History_Population.upper_decs,Elite_Population.upper_decs),[],2);
            ID1 = find(Dis_temp>0);
            ID2 = find(Dis_temp==0);
            ID = Candidate_ID(History_Population(ID1),Elite_Population);
            
            Elite_Population_temp = History_Population([ID1(ID);ID2]);
            
            [~,History_data] = lower_level_optimize(Elite_Population_temp,History_data,redundant_dec_lobj,T,alpha,belta,Maxgen);
            
            [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,redundant_dec_lobj,W,T,alpha,belta,Maxgen,1);

            Elite_Population_temp = History_data.bilevel_best;
            Dis = min(pdist2(Elite_Population.upper_decs,Elite_Population_temp.upper_decs),[],2);
            Elite_Population(Dis==0) = [];
            Elite_Population = cat(2,Elite_Population,Elite_Population_temp);
    end
    
    if length(Elite_Population) > N
        Elite_Population = Upper_Select(Elite_Population,N,W,'bilevel');
    end
    
    
end

function ID = Candidate_ID(Elite_Population_temp,Elite_Population)
    ID = [];
    if isempty(Elite_Population_temp)
        return
    end
    
    Dis = min(pdist2(Elite_Population_temp.upper_decs,Elite_Population.upper_decs),[],2);
    
    ind1 = find(Dis>0);
    if ~isempty(ind1)
        
        [I_epsilon,rate] = C_I_epsilon(Elite_Population_temp(ind1).upper_objs,Elite_Population.upper_objs);
        ID = ind1(I_epsilon>=1|rate>=0.8*length(ind1));
        
        ID = unique(ID);
    end
end

function [I_epsilon,rate] = C_I_epsilon(PopObj1,PopObj2)

    N1 = size(PopObj1,1);
    N2 = size(PopObj2,1);
    Min_obj = min([PopObj1;PopObj2],[],1);
    
    PopObj1 = PopObj1 - Min_obj;
    PopObj2 = PopObj2 - Min_obj;
    
    I_epsilon = zeros(1,N1);
    rate = zeros(1,N1);
    for i=1:N1
        I_epsilon(i) = min(max(PopObj2./max(repmat(PopObj1(i,:),N2,1),eps),[],2));
        rate(i) = sum(max(repmat(PopObj1(i,:),N2,1)./max(PopObj2,eps),[],2)<1);
    end
    
end