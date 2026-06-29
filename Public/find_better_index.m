function I = find_better_index(Pop,individual,type)
    Pop_temp = [Pop,individual];
    switch type
        case 'upper'
            FrontN = NDSort(Pop_temp.upper_objs,Pop_temp.upper_CVs,inf);
            I = find(FrontN(1:end-1)<FrontN(end));
        case 'lower'
            FrontN = NDSort(Pop_temp.lower_objs,Pop_temp.lower_CVs,inf);
            I = find(FrontN(1:end-1)<FrontN(end));
    end
end