function Upper_Data_Population = Upper_Data(Pop)

    Upper_Data_Population = struct('ID',[],'Population',[],'FrontNo',[],'CrowdDis',[]);
    Population = [];
    ID = [];
    for i=1:length(Pop)
        Population_i = Pop(i).Population;
        Population = cat(2,Population,Population_i);
        ID = cat(2,ID,repmat(i,1,length(Population_i)));
    end
    PopObj = Population.upper_objs;
    FrontNo = NDSort(PopObj,Population.upper_cons,inf);
    CrowdDis = CrowdingDistance(PopObj,FrontNo);
    
    Upper_Data_Population.ID = ID;
    Upper_Data_Population.Population = Population;
    Upper_Data_Population.FrontNo = FrontNo;
    Upper_Data_Population.CrowdDis = CrowdDis;
end