function History_data = Update_Data(History_data,Elite_Population,N)
    
    Elite_upper_decs = unique(Elite_Population.upper_decs,'row');
    [Dis,inds] = min(pdist2(Elite_upper_decs,History_data.upper_decs),[],2);
    Elite_data_ID = inds(Dis==0); 
    
    LH = length(History_data);
    if LH > N
        [History_Population,ID] = History_data.Best_Archives;
        FID = History_Population.bilevel_feasibles;
        FPop = History_Population(FID);
        
        ID = ID(FID);
        UID = unique(ID);
        L_UID = length(UID);
        
        if  L_UID >= N
            [FrontNo,MaxNo] = NDSort(FPop.upper_objs,inf);
            i = 1;
            Next = false(1,LH);
            while sum(Next)<N && i<= MaxNo
                ID_i = unique(ID((FrontNo==i)));
                if sum(Next)+sum(Next(ID_i)==false)<=N
                    Next(ID_i) = true;
                else
                    ID_i = ID_i(Next(ID_i)==false);
                    Next(ID_i(randperm(length(ID_i),N-sum(Next)))) = true;
                end
                i=i+1;
            end
            Next(Elite_data_ID) = true;
            History_data = History_data(Next);
            
        else
            Data = History_data(UID);
            P = 1:LH;
            P(UID) = [];
            [inFPop,I] = History_data(P).Best_Archives;
            [~,rank] = sort(inFPop.upper_CVs,'ascend');
            History_data = cat(2,Data,History_data(P(I(rank(1:N-L_UID)))));
        end
    end
    
end