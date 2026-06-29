function [History_data,History_Population] = Update_Data(History_data,N)
    
    LH = length(History_data);
    Next =false(1,LH);
    History_Population = [];
    
    [Population,ID] = History_data.Best_Archives;
    [FrontNo,MaxNo] = NDSort(Population.upper_objs,Population.cons,inf);
    CrowdDis = CrowdingDistance(Population.upper_objs,FrontNo);
    for i = 1:MaxNo
        if sum(Next)<N
            ID_i = unique(ID((FrontNo==i)));
            ID_i = ID_i(Next(ID_i)==false);
            if sum(Next)+length(ID_i)<=N
                Next(ID_i) = true;
            else
                CrowdDis_i = [];
                for j = ID_i
                    index = find(FrontNo == i & ID==j);
                    [~,Rank_j] = sort(CrowdDis(FrontNo == i & ID==j ),'descend');
                    CrowdDis_i = cat(2,CrowdDis_i,CrowdDis(index(Rank_j(1))));
                end
                [~,Rank] = sort(CrowdDis_i,'descend');
                ID_i = ID_i(Rank(1:N-sum(Next)));
                Next(ID_i) = true;
            end
            for j = ID_i
                index = find(FrontNo == i & ID==j);
                [~,Rank_j] = sort(CrowdDis(index),'descend');
                History_Population = cat(2,History_Population,Population(index(Rank_j(1))));
            end
        else
            break;
        end
    end
    
    History_data = History_data(Next);

end