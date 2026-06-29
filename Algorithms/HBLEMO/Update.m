function [Pop_Next,Upper_Data_Pop_Next,ID] = Update(Pop,Off,N_u)
    N=length(Pop);
    Pop = cat(2,Pop,Off);
    Upper_Data_Pop = Upper_Data(Pop);
    
    ID_Next = [];
    FrontNo = Upper_Data_Pop.FrontNo;
    ID = Upper_Data_Pop.ID;
    
    for i=1:max(FrontNo)
        if length([Pop(ID_Next).Population])<N_u
            Pop_i = Upper_Data_Pop.Population(FrontNo == i);
            ID_i = ID(FrontNo == i);
            unique_ID_i = unique(ID_i);
            for j = unique_ID_i
                if ~ismember(j,ID_Next)
                    Dis =  min(pdist2(Pop_i(ID_i==j).lower_decs,Pop(j).Population(Pop(j).FrontNo ==1).lower_decs),[],2);
                    if any(Dis==0)
                        ID_Next = unique(cat(2,ID_Next,j));
                    end
                end
            end
        else
            break;
        end
    end
    
    Pop_Next = Pop(ID_Next);
    ID = find(ID_Next<=N);
    
    Upper_Data_Pop_Next = Upper_Data(Pop_Next);
end