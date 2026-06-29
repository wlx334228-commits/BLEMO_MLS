function Offspring = DE_current_2_bin(Target,Parent1,Parent2,Parent3,Parent4,type,Parameter)
%DE - Differential evolution operator.
    Global=GLOBAL.GetObj();
    if isempty(Target)
        Offspring = [];
        return
    end
    %% Parameter setting
    if nargin > 6
        [CR,F,proM,disM] = deal(Parameter{:});
    else
        [CR,F,proM,disM] = deal(1,0.5,1,20);
    end
    if isa(Target,'INDIVIDUAL')
        calObj = true;
        switch type
            case 'upper'
                Target = Target.upper_decs;
                Parent1 = Parent1.upper_decs;
                Parent2 = Parent2.upper_decs;
                Parent3 = Parent3.upper_decs;
                Parent4 = Parent4.upper_decs;
            case 'lower'
                Target = Target.lower_decs;
                Parent1 = Parent1.lower_decs;
                Parent2 = Parent2.lower_decs;
                Parent3 = Parent3.lower_decs;
                Parent4 = Parent4.lower_decs;
        end
    else
        calObj = false;
    end
    
    %% DE Mutation and Crossover
    [N,D]  = size(Target);
    Site = rand(N,D) < CR;
    K = randi(D,N,1);
    for i=1:N
        Site(i,K(i)) = true;
    end
    Offspring = Target;
    Offspring(Site) = Target(Site) + rand*F*(Parent2(Site) - Parent1(Site)) +  F*(Parent4(Site) - Parent3(Site));
    
    Offspring = Global.problem.Decs(Offspring,type);
    
    %% Polynomial mutation
%     if proM>0
%         switch type
%             case 'upper'
%                 Lower = repmat(Global.upper_domain(1,:),N,1);
%                 Upper = repmat(Global.upper_domain(2,:),N,1);
%             case 'lower'
%                 Lower = repmat(Global.lower_domain(1,:),N,1);
%                 Upper = repmat(Global.lower_domain(2,:),N,1);
%             case 'upper&lower'
%                 Lower = repmat([Global.upper_domain(1,:),Global.lower_domain(1,:)],N,1);
%                 Upper = repmat([Global.upper_domain(2,:),Global.lower_domain(2,:)],N,1);
%         end
%         Site  = rand(N,D) < proM/D;
%         mu    = rand(N,D);
%         temp  = Site & mu<=0.5;
%         Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%             (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%         temp = Site & mu>0.5;
%         Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%             (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%         Offspring = Global.problem.Decs(Offspring,type);
%     end
    
%     if calObj
%         Offspring = INDIVIDUAL(Offspring,type);
%     end
end
