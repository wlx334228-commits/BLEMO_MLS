function BLPIEA(Global)
    %% Parameter setting
    % K --- the zone number for M2M deviding;
    % T --- the size of the neighbourhood;
    % alpha --- the period for group lower search;
    % beta --- a parameter for termination check in group lower search;
    % theta -- a parameter for termination check in upper search;
    % UMaxG --- the maximum generation for detective upper search;
    % LMaxG1 --- the maximum generation for group lower search;
    % LMaxG2 --- the maximum generation for local group lower search;
    
    [K,T1,alpha,beta,theta,LMaxG] = Global.ParameterSet(10,10,10,10,1e-3,500);
    History_data = [];
    Elite_Population = [];
    
    %% Generate the reference points and random population
    upper_N = Global.N(1);
    W = UniformPoint(K,Global.M(1));
    N = Global.N(3);
    
    %% Generate random population
    Population = Global.Initialization;
    Population = Global.Evaluate(Population,'upper');
    
    [Population,History_data] = lower_level_optimize(Population,History_data,T1,beta,LMaxG);
    
    Explored_Population = tentative_explore(Population,beta,theta);
    [~,History_data] = lower_level_optimize(Explored_Population,History_data,T1,beta,LMaxG);
     
    [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,W,T1,beta,LMaxG,1);
    
    Elite_Population_temp = Elite_Population;
    
    Population = Upper_Select(History_data.Best_Archives,upper_N,W,'bilevel');
    Global.upper_Output(Population);
    
    while Global.NotTermination(Elite_Population,true)
        %% Global upper search combined with group lower search
        labels = Population.labels;
        ID = find(labels);
        [Population(ID),History_data] = lower_level_optimize(Population(ID),History_data,T1,beta,LMaxG);
        History_Population = Upper_Select(History_data.Best_Archives,upper_N,W,'bilevel');
        [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,W,T1,beta,LMaxG,1);
        
        NDrate = length(History_Population.upper_best)/length(History_Population);
        
        Distance = sort(pdist2(History_Population.upper_objs,History_Population.upper_objs),2);
        D = max(Distance(:,1:1+ceil(0.5*(1-NDrate)*length(History_Population))),[],2);
        
%         idx = kmeans(History_Population.upper_objs,10)
        
        Offspring = Generate(Population,History_Population,D,W,T1);
        
        if mod(Global.gen,alpha) == 0
            
            [~,History_data] = lower_level_optimize([Population,Offspring],History_data,T1,beta,LMaxG);
            Population = Upper_Select(History_data.Best_Archives,upper_N,W,'bilevel');
            
            Elite_Population = Elite_Population_temp;
            History_data = Update_Data(History_data,Elite_Population,N);
            
            Elite_Population_temp = History_data.bilevel_best;
            if length(Elite_Population_temp) > N
                Elite_Population_temp = Upper_Select(Elite_Population_temp,N,W,'bilevel');
            end
            
        else
            labels = Population.labels;
            Offspring= [Population(~labels),Offspring];
            [Distance,ID] = min(pdist2(Offspring.upper_objs,History_Population.upper_objs),[],2);
            I = sum(Distance>D(ID))/length(Distance);
            if rand < I
                if sum(Distance<=D(ID))>=upper_N-sum(labels)
                    Offspring = Offspring(Distance<=D(ID));
                else
                    inds = find(Distance>D(ID));
                    [~,rank] = sort(abs(Distance(inds)-D(ID(inds))),'ascend');
                    Offspring = [Offspring(Distance<=D(ID)),Offspring(inds(rank(1:upper_N-(sum(labels)+sum(Distance<=D(ID))))))];
                end
            end
            
            Population = Upper_Select([Population(labels),Offspring],upper_N,W,'upper');
            
            [Elite_Population_temp,History_data] = Update_Elites(Elite_Population_temp,History_data,W,T1,beta,LMaxG,2);
            [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,W,T1,beta,LMaxG,1);
        end
        
        if mod(Global.gen+1,alpha) == 0
            if any([Global.upper_FEs,Global.lower_FEs] > Global.maxFEs)
                Global.NotTermination(Elite_Population,false);
            end
        end
                
        Global.upper_Output(Population,2);
        
    end
end
