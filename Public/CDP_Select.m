function Population = CDP_Select(Population,N,type)
    switch type
        case 'upper'
            Pop_objs = Population.upper_objs;
            Pop_cons = Population.upper_cons;
        case 'lower'
            Pop_objs = Population.lower_objs;
            Pop_cons = Population.lower_cons;
    end
    
    CV = sum(max(0,Pop_cons),2);
    F_index = find(CV==0);
    if length(F_index)>N
        [~,rank]=sort(Pop_objs(F_index),'ascend');
        Population = Population(F_index(rank(1:N)));
    else
        [~,rank]=sort(CV,'ascend');
        Population = Population(rank(1:N));
    end
end