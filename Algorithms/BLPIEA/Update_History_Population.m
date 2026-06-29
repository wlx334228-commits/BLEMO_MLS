function History_Population = Update_History_Population(History_Population,History_data)
    RV = History_Population.lower_objs;
    RV = RV./sqrt(sum(RV.^2,2));
    
    [unique_upper_decs,IA] = uniquetol(History_Population.upper_decs,0,'ByRows', true,'OutputAllIndices', true);
    [Dis,Ind] = min(pdist2(unique_upper_decs,History_data.upper_decs),[],2);
    index = find(Dis==0);
    History_Population = Update(History_Population,History_data(Ind(index)),IA(index),RV);
end

function Population = Update(Population,Data,IA,RV)
    Global = GLOBAL.GetObj();
    S = length(Data);
    for i = 1:S
        Best_Archive_i = Data(i).Best_Archive;
        L = length(Best_Archive_i);
        ID = IA{i};
        K = length(ID);
        if L>K
            Next = false(1,L);
            for j=1:K
                index = find(~Next);
                [~,pos] = min(max(Best_Archive_i(index).lower_objs./RV(ID(j),:),[],2));
                Next(index(pos)) = true;
            end
            Population(ID) = Best_Archive_i(Next);
%             [~,inds] = min(pdist2(Population(ID).lower_decs,Best_Archive_i.lower_decs),[],2);
%             Population(ID) = Best_Archive_i(inds);
%             ID = ID(dis~=0);
%             if ~isempty(ID)
%                 K = length(ID);
%                 Best_Archive_i(inds(dis==0)) = [];
%                 L = length(Best_Archive_i);
%                 if L <= K
%                     Population(ID) = [Best_Archive_i,Best_Archive_i(randi(L,1,K-L))];
%                 else
%                     for j=1:K
%                         [~,inds] = min(pdist2(Population(ID(j)).lower_decs,Best_Archive_i.lower_decs),[],2);
%                         Population(ID(j)) = Best_Archive_i(inds);
%                         Best_Archive_i(inds) = [];
%                     end
%                 end
%             end
        else
            Population(ID) = [Best_Archive_i,Best_Archive_i(randi(L,1,K-L))];
        end
    end
    
end