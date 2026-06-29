function [Population,Next] = Upper_Select(Population,N,W,type)
    if length(Population)>N
        Global = GLOBAL.GetObj();
        switch type
            case 'upper'
                Population = Global.Evaluate(Population,'upper');
                feasible_index = find(Population.upper_feasibles);
            case 'lower'
                Population = Global.Evaluate(Population,'lower');
                feasible_index = find(Population.lower_feasibles);
            case 'bilevel'
                
                Population = Global.Evaluate(Population,'upper');
                feasible_index = find(Population.bilevel_feasibles);
        end
        
        if length(feasible_index)>N
            FPop = Population(feasible_index);
            popsize = length(FPop);
            K = size(W,1);
            S = N/K;
            
            B = 1-pdist2(W,W,'cosine');
            [~,B] = sort(B,2,'descend');
            
            switch type
                case 'upper'
                    PopObj = FPop.upper_objs;
                case 'lower'
                    PopObj = FPop.lower_objs;
                case 'bilevel'
                    PopObj = FPop.upper_objs;
            end
            
            Zmin = min(PopObj,[],1);
            Zmax = max(PopObj,[],1);
            PopObj = (PopObj-repmat(min(PopObj,[],1),popsize,1))./(repmat(max(Zmax-Zmin,1e-6),popsize,1));
            
            [~,ZoneNum] = max(1-pdist2(PopObj,W,'cosine'),[],2);
            
            Next = [];
            current = [];
            for i = 1:K
                current = cat(1,current,ZoneNum'==i);
            end
            current_N = sum(current,2);
            
            for i=1:K
                t=1;
                while sum(current_N(B(i,1:t)))<S
                    t = t + 1;
                end
                currentID = find(ZoneNum'==i);
                neighborID = [];
                for j = B(i,2:t)
                    neighborID = cat(2,neighborID,find(ZoneNum'==j));
                end
                if length(currentID) == S
                    Next = cat(2,Next,currentID);
                elseif length(currentID) > S
                    [~,IndexP] = NSGAII_Update(FPop(currentID),S,'upper');
                    Next = cat(2,Next,currentID(IndexP));
                elseif length(currentID)+length(neighborID)==S
                    Next = cat(2,Next,[currentID,neighborID]);
                else
                    [~,IndexP] = NSGAII_Update(FPop(neighborID),S-length(currentID),'upper');
                    Next = cat(2,Next,[currentID,neighborID(IndexP)]);
                end
                ZoneNum(Next) = 0;
                current(:,Next) = 0;
                current_N = sum(current,2);
            end
            
            Next = feasible_index(Next);
            Population = Population(Next);
        else
            Next = false(1,length(Population));
            Next(feasible_index) = true;
            Inds = find(~Next);
            switch type
                case 'bilevel'
                    CV = Population(Inds).upper_CVs;
                case 'upper'
                    CV = Population(Inds).upper_CVs;
                case 'lower'
                    CV = Population(Inds).lower_CVs;
            end
            [~,rank] = sort(CV,'ascend');
            Next(Inds(rank(1:N-length(feasible_index))))=true;
            Population = Population(Next);
        end
    else
        Next = 1:length(Population);
    end
   
end