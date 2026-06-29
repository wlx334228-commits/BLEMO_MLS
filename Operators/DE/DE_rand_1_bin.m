function Offspring = DE_rand_1_bin(Target,Parent1,Parent2,Parent3,type,Parameter)
%DE - Differential evolution operator.
    Global=GLOBAL.GetObj();
    if isempty(Target)
        Offspring = [];
        return
    else
        Offspring = Target;
    end
    %% Parameter setting
    if nargin > 4
        [CR,F] = deal(Parameter{:});
    else
        [CR,F] = deal(1,0.5,1,20);
    end
    if isa(Parent1(1),'INDIVIDUAL')
        calObj = true;
        switch type
            case 'upper'
                Target = Target.upper_decs;
                Parent1 = Parent1.upper_decs;
                Parent2 = Parent2.upper_decs;
                Parent3 = Parent3.upper_decs;
            case 'lower'
                Target = Target.upper_decs;
                Parent1 = Parent1.lower_decs;
                Parent2 = Parent2.lower_decs;
                Parent3 = Parent3.lower_decs;
        end
    else
        calObj = false;
    end
    
    %Mutation
    [N,D]  = size(Parent1);
    trial = Parent1 + F*(Parent2 - Parent3);
    
    %Crossover
    Site = rand(N,D) < CR;
    K = randi(D,N,1);
    for i=1:N
        Site(i,K(i)) = true;
    end
    
    Offspring = Target;
    Offspring(Site) = trial(Site);
    
    Offspring = Global.problem.Decs(Offspring,type);
    
%     if calObj
%         Offspring = INDIVIDUAL(Offspring,type);
%     end
end
