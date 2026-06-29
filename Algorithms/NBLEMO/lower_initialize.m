function Population = lower_initialize(upper_dec,Upper_Data_Pop,Archive,N_l_0,N_l)
    Global = GLOBAL.GetObj();
    
    Population = [];
    Data_i.upper_dec = upper_dec;
    lower_decs = [];
    
    if N_l==N_l_0
        for j=1:N_l
            if rand < length(Archive.Population)/(length(Archive.Population)+length(Upper_Data_Pop.Population))
                MatingPool = TournamentSelection(2,2,-Archive.CrowdDis);
                lower_dec_j = GAhalf(Archive.Population(MatingPool),'lower',{0.9,15,0.1,20});
            else
                MatingPool = TournamentSelection(2,2,Upper_Data_Pop.FrontNo,-Upper_Data_Pop.CrowdDis);
                lower_dec_j = GAhalf(Upper_Data_Pop.Population(MatingPool),'lower',{0.9,15,0.1,20});
            end
            lower_decs = cat(1,lower_decs,lower_dec_j);
        end
    else
        [~,ID] = min(pdist2(upper_dec,Archive.Population.upper_decs),[],2);
        Dis = min(pdist2(Upper_Data_Pop.Population.upper_decs,Archive.Population(ID).upper_decs),[],2);
        lower_decs = Upper_Data_Pop.Population(Dis==0).lower_decs;
        if size(lower_decs,1)>N_l
            lower_decs = lower_decs(randperm(size(lower_decs,1),N_l),:);
        else
            for j=1:N_l-size(lower_decs,1)
                if rand < length(Archive.Population)/(length(Archive.Population)+length(Upper_Data_Pop))
                    MatingPool = TournamentSelection(2,2,-Archive.CrowdDis);
                    lower_dec_j = GAhalf(Archive.Population(MatingPool),'lower',{0.9,15,0.1,20});
                else
                    MatingPool = TournamentSelection(2,2,Upper_Data_Pop.FrontNo,-Upper_Data_Pop.CrowdDis);
                    lower_dec_j = GAhalf(Upper_Data_Pop.Population(MatingPool),'lower',{0.9,15,0.1,20});
                end
                lower_decs = cat(1,lower_decs,lower_dec_j);
            end
        end
    end
    
    Data_i.Population = Global.Evaluate(INDIVIDUAL({upper_dec,lower_decs}),'lower');
    [~,~,Data_i.FrontNo,Data_i.CrowdDis] = NSGAII_Update(Data_i.Population,N_l,'lower');
    Population = cat(2,Population,Data_i);
end