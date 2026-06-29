function [Population,IndexP,FrontNo,CrowdDis] = NSGAII_Update(Population,N,type)
    switch type
        case 'bilevel'
            PopObj = Population.upper_objs;
            [FrontNo,MaxFNo] = NDSort(PopObj,[Population.upper_cons,Population.lower_cons],N);
            CrowdDis = CrowdingDistance(PopObj,FrontNo);
        case 'upper'
            PopObj = Population.upper_objs;
            [FrontNo,MaxFNo] = NDSort(PopObj,Population.upper_cons,N);
            CrowdDis = CrowdingDistance(PopObj,FrontNo);
        case 'lower'
            PopObj = Population.lower_objs;
            [FrontNo,MaxFNo] = NDSort(PopObj,Population.lower_cons,N);
            CrowdDis = CrowdingDistance(PopObj,FrontNo);
    end
    
    IndexP = FrontNo < MaxFNo;
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    IndexP(Last(Rank(1:N-sum(IndexP)))) = true;
    
    %% Population for next generation
    Population = Population(IndexP);
    IndexP = find(IndexP);
    FrontNo = FrontNo(IndexP);
    CrowdDis = CrowdDis(IndexP);
end