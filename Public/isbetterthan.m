function I = isbetterthan(Pop1,Pop2,type)
    Pop = [Pop1,Pop2];
    switch type
        case 'upper'
            FrontN = NDSort(Pop.upper_objs,Pop.upper_CVs,1);
            if FrontN(1) < FrontN(2)
                I = true;
            else
                I = false;
            end
        case 'lower'
            FrontN = NDSort(Pop.lower_objs,Pop.lower_CVs,1);
            if FrontN(1) < FrontN(2)
                I = true;
            else
                I = false;
            end
    end
end