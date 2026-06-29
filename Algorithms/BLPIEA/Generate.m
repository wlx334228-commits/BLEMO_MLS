function Offspring = Generate(Population,History_Population,D,W,T)
    Global = GLOBAL.GetObj();
    N = length(Population);
    
    [~,New_Individuals] = Population.labels(false);
    if ~isempty(New_Individuals)
        [Distance,ID] = min(pdist2(New_Individuals.upper_objs,History_Population.upper_objs),[],2);
        I = sum(Distance>D(ID))/length(Distance);
%         MatingPop = NSGAII_Update([History_Population,New_Individuals(Distance<=D(ID))],N,'upper');
        MatingPop = [History_Population,New_Individuals(Distance<=D(ID))];
    else
        MatingPop = Population;
        I = 0;
    end
    
    [Target,Pool] = Upper_Mating(MatingPop,N,W,T,2);
%     if rand < I
%         upper_decs = DE_current_1_bin(MatingPop(Target).upper_decs,MatingPop(Pool(:,1)).upper_decs,...
%                      MatingPop(Pool(:,2)).upper_decs,'upper',{0.5,0.2,0.1,20});
%     else
        upper_decs = DE_current_1_bin(MatingPop(Target).upper_decs,MatingPop(Pool(:,1)).upper_decs,...
                     MatingPop(Pool(:,2)).upper_decs,'upper',{0.5,0.5,0.1,20});
%     end
    [~,LinkID] = min(pdist2(upper_decs,History_Population.upper_decs),[],2);
    Target = Global.Evaluate(INDIVIDUAL({upper_decs,History_Population(LinkID).lower_decs}),'upper');
    
%     [~,LinkID] = min(pdist2(LinkPop.upper_objs,MatingPop.upper_objs),[],2);
%     MatingPop = cat(2,History_Population,MatingPop(unique(LinkID)));
    
%     Pool = [];
%     for i=1:N
% %         if rand <I
% %             Pool = cat(1,Pool,History_Population(randperm(length(History_Population),2)));
% %         else
%             Pool = cat(1,Pool,MatingPop(randperm(length(MatingPop),2)));
% %         end
%     end
%     lower_decs = DE_current_1_bin(LinkPop.lower_decs,Pool(:,1).lower_decs,...
%                  Pool(:,2).lower_decs,'lower',{0.5,0.2,0.1,20});    

    B = pdist2(Target.upper_objs,MatingPop.upper_objs);
    B(B==0) = inf;
    [~,B] = sort(B,2);
    B = B(:,1:T);
    Pool = [];
    for i=1:N
        if rand < 0.5
            Pool = cat(1,Pool,B(i,randperm(T,2)));
        else
            Pool = cat(1,Pool,randperm(N,2));
        end
    end
    
    lower_decs = DE_current_1_bin(Target.lower_decs,MatingPop(Pool(:,1)).lower_decs,...
                 MatingPop(Pool(:,2)).lower_decs,'lower',{0.5,0.2,0,20});
    
    Offspring = INDIVIDUAL({upper_decs,lower_decs});
    Offspring = Global.Evaluate(Offspring,'upper');
end
