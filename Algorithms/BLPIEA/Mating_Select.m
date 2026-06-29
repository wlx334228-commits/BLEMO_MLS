function MatingPop = Mating_Select(Population,History_Population,N,W)
    Global = GLOBAL.GetObj();
    [~,labels] = Population.labels;
    L = length(Population);
    [FrontNo,MaxFNo] = NDSort(Population.upper_objs,Population.upper_cons,inf);
    if any(labels)
        MaxFNo = min(FrontNo(labels));
        L = sum(FrontNo<=MaxFNo);
    end
    for i = MaxFNo:-1:1
        if sum(FrontNo>=i & FrontNo<=MaxFNo) >= L/2
            break;
        end
    end
    
    Population = [History_Population,Population(FrontNo>=i)];
    MatingPop = NSGAIII_Update(Population,N,'upper');
end
    

    