function HBLEMO(Global)
    tao = 10;
    Epsilon_u = 1e-2;
    Epsilon_l = 0.1;
    N_u = Global.N(3);%20*sum(Global.D);
    N_l_0 = Global.N(2);%round(sqrt(N_u*Global.D(2)/Global.D(1)));
    n_s_0 = Global.N(1);%round(sqrt(N_u*Global.D(1)/Global.D(2)));
    N_A = 10*N_u;
    
    Archive = struct('Population',[],'CrowdDis',[]);
    
    [Population,Upper_Data_Pop] = global_initialize(n_s_0,N_l_0);
    
    [Population,Archive,~,T_l_max] = ll_NSGAII(Population,Upper_Data_Pop,Archive,N_l_0,N_A,tao,Epsilon_l);
    
    Archive_Population = Archive.Population;
    
    Upper_Data_Pop = Upper_Data(Population);
    
    Global.upper_Output(Upper_Data_Pop.Population,2);
    
    Objs = {};
    
    while Global.NotTermination(Archive_Population,true)
        
        N_Q = 0;
        Offspring =[];
        
        while N_Q < N_u
            
            if rand < length(Archive.Population)/(length(Archive.Population)+length(Upper_Data_Pop.Population))
                MatingPool = TournamentSelection(2,2,-Archive.CrowdDis);
                upper_dec = GAhalf(Archive.Population(MatingPool),'upper',{0.9,15,0.1,20});
            else
                MatingPool = TournamentSelection(2,2,Upper_Data_Pop.FrontNo,-Upper_Data_Pop.CrowdDis);
                upper_dec = GAhalf(Upper_Data_Pop.Population(MatingPool),'upper',{0.9,15,0.1,20});
            end
            
            if ~isempty(Archive.Population)
                Archive_upper_decs = Archive.Population.upper_decs;
                Delta_U = max(max(pdist2(Archive_upper_decs,Archive_upper_decs),[],'all'),eps);
                Delta_u = min(pdist2(upper_dec,Archive_upper_decs),[],2);
                N_l = min(max(round(N_l_0*Delta_u/Delta_U),4),N_l_0);
                T_l = min(T_l_max*Delta_u/Delta_U,T_l_max);
            else
                N_l = N_l_0;
                T_l = T_l_max;
            end
            
            individual = lower_initialize(upper_dec,Upper_Data_Pop,Archive,N_l_0,N_l);
            [individual,Archive,N_individual] = ll_NSGAII(individual,Upper_Data_Pop,Archive,N_l_0,N_A,tao,Epsilon_l,T_l);
            
            N_Q = N_Q + N_individual;
            
            Offspring = cat(2,Offspring,individual);
            
        end
        
        
        [Population,Upper_Data_Pop,ID] = Update(Population,Offspring,N_u);
        
        [Population(ID),Archive] = ll_NSGAII(Population(ID),Upper_Data_Pop,Archive,N_l_0,N_A,tao,Epsilon_l,T_l_max);
        
        Upper_Data_Pop = Upper_Data(Population);
        
        Global.upper_Output(Upper_Data_Pop.Population,2);
        
%         Objs = cat(2,Objs,{Upper_Data_Pop.Population(Upper_Data_Pop.FrontNo==1).upper_objs});
        if length(Archive.Population)>N_u
            [~,rank] = sort(Archive.CrowdDis,'descend');
            Archive_Population = Archive.Population(rank(1:N_u));
        else
            Archive_Population = Archive.Population;
        end
        
        if all([Global.upper_FEs,Global.lower_FEs] > Global.maxFEs)
            Global.NotTermination(Archive.Population,false);
        end
        
%         if ~isempty(Archive_Population)    
%             Objs = cat(2,Objs,{Archive_Population.upper_objs});
%         end
%         
%         if length(Objs)==tao
%             RefPoint = max(cat(1,Objs{:}),[],1);
%             HV = cell2mat(cellfun(@(x) CalHV(x,RefPoint), Objs, 'UniformOutput',false));
%             H_u = (max(HV)-min(HV))/(max(HV)+min(HV));
%             if H_u <= Epsilon_u && all([Global.upper_FEs,Global.lower_FEs] > Global.maxFEs)
%                 Global.NotTermination(Archive_Population,false);
%             end
%             Objs(1) = [];
%         end
    end
end