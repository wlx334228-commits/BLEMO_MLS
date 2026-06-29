function Offspring = DE_rand_to_best_1(Parent1,Parent2,Parent3,best_dec,type,Parameter)
%DE - Differential evolution operator.
    if isempty(Parent1)
        Offspring = [];
        return
    end
    %% Parameter setting
    if nargin > 5
        [CR,F1,F2,proM,disM] = deal(Parameter{:});
    else
        [CR,F1,F2,proM,disM] = deal(1,0.5,0.5,1,20);
    end
    if isa(Parent1,'INDIVIDUAL')
        calObj = true;
        switch type
            case 'upper'
                Parent1 = Parent1.upper_decs;
                Parent2 = Parent2.upper_decs;
                Parent3 = Parent3.upper_decs;
                best_dec = best_dec.upper_decs;
            case 'lower'
                Parent1 = Parent1.lower_decs;
                Parent2 = Parent2.lower_decs;
                Parent3 = Parent3.lower_decs;
                best_dec = best_dec.lower_decs;
        end
    else
        calObj = false;
    end
    
    [N,D]  = size(Parent1);
    if size(best_dec,1)==1
        best_dec = repmat(best_dec,N,1);
    end
    Site = rand(N,D) < CR;
    Offspring = Parent1;
    Offspring(Site) = Offspring(Site) + F1*(best_dec(Site)-Parent1(Site)) + F2 *(Parent2(Site)-Parent3(Site));
    
    %% Polynomial mutation
    Global=GLOBAL.GetObj();
    switch type
        case 'upper'
            Lower = repmat(Global.upper_domain(1,:),N,1);
            Upper = repmat(Global.upper_domain(2,:),N,1);
        case 'lower'
            Lower = repmat(Global.lower_domain(1,:),N,1);
            Upper = repmat(Global.lower_domain(2,:),N,1);
        case 'upper&lower'
            Lower = repmat([Global.upper_domain(1,:),Global.lower_domain(1,:)],N,1);
            Upper = repmat([Global.upper_domain(2,:),Global.lower_domain(2,:)],N,1);
    end
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%     if calObj
%         Offspring = INDIVIDUAL(Offspring,type);
%     end
end