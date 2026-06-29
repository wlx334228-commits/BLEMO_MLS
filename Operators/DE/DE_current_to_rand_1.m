function Offspring = DE_current_to_rand_1(Target,Parent1,Parent2,Parent3,type,Parameter)
%DE - Differential evolution operator.
    Global=GLOBAL.GetObj();
    %% Parameter setting
    if nargin > 5
        [CR,F] = deal(Parameter{:});
    else
        [CR,F] = deal(1,0.5);
    end
    
    if isempty(Target)
        Offspring = [];
        return
    else
        Offspring = Target;
    end
    
    if isa(Target,'INDIVIDUAL')
        calObj = true;
        switch type
            case 'upper'
                Target = Target.upper_decs;
                Parent1 = Parent1.upper_decs;
                Parent2 = Parent2.upper_decs;
                Parent3 = Parent3.upper_decs;
            case 'lower'
                Target = Target.lower_decs;
                Parent1 = Parent1.lower_decs;
                Parent2 = Parent2.lower_decs;
                Parent3 = Parent3.lower_decs;
        end
    else
        calObj = false;
    end
    
    Offspring = Target + rand*(Parent1-Target) + rand*F*(Parent2-Parent3);
    
    Offspring = Global.problem.Decs(Offspring,type);

    
%     if calObj
%         Offspring = INDIVIDUAL(Offspring,type);
%     end
end