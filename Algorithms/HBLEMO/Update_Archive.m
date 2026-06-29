function Archive = Update_Archive(Archive,Archive_P,N_A)
    Archive.Population = cat(2,Archive.Population,Archive_P.Population);
    PopObj = Archive.Population.upper_objs;
    [FrontNo,~] = NDSort(PopObj,Archive.Population.upper_cons,1);
    if sum(FrontNo==1) <= N_A
        CrowdDis = CrowdingDistance(PopObj,FrontNo);
        Archive.Population = Archive.Population(FrontNo==1);
        Archive.CrowdDis = CrowdDis(FrontNo==1);
    else
        [Archive.Population,~,~,Archive.CrowdDis] = NSGAII_Update(Archive.Population(FrontNo==1),N_A,'lower');
    end
end