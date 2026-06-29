function [best,index] = find_best(Pop,level)
    switch level
        case 'upper'
%             FrontNo = NDSort(Pop.upper_objs,Pop.upper_CVs,1);
%             best = Pop(FrontNo == 1);
%             index = find(FrontNo == 1);
            feasible = Pop.upper_feasibles;
            CV = Pop.upper_CVs;
            F_index = find(feasible==1);
            FSet = Pop(F_index);
            if ~isempty(FSet) > 0
                [~,ind] = min(FSet.upper_objs);
                best = FSet(ind);
                index = F_index(ind);
            else
                [~,ind] = min(CV);
                best = Pop(ind);
                index = ind;
            end
        case 'lower'
%             FrontNo = NDSort(Pop.lower_objs,Pop.lower_CVs,1);
%             best = Pop(FrontNo == 1);
%             index = find(FrontNo == 1);
            feasible = Pop.lower_feasibles;
            CV = Pop.lower_CVs;
            F_index = find(feasible==1);
            FSet = Pop(F_index);
            if ~isempty(FSet) > 0
                [~,ind] = min(FSet.lower_objs);
                best = FSet(ind);
                index = F_index(ind);
            else
                [~,ind] = min(CV);
                best = Pop(ind);
                index = ind;
            end
    end
end