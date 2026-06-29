function Offspring = Polynomial_mutation(Parent,Lower,Upper,type,Parameter)
% PCX

    %% Parameter setting
    if nargin > 2
        [proM,disM] = deal(Parameter{:});
    else
        [proM,disM] = deal(1,20);
    end
    if isa(Parent,'INDIVIDUAL')
%         calObj = true;
        switch type
            case 'upper'
                Parent = Parent.upper_decs;
            case 'lower'
                Parent = Parent.lower_decs;
        end
%     else
%         calObj = false;
    end
    
    [N,D]  = size(Parent);
    Offspring = Parent;
    
    %% Polynomial mutation
    Lower = repmat(Lower,N,1);
    Upper = repmat(Upper,N,1);
    
%     Global=GLOBAL.GetObj();
%     switch type
%         case 'upper'
%             Lower = repmat(Global.upper_domain(1,:),N,1);
%             Upper = repmat(Global.upper_domain(2,:),N,1);
%         case 'lower'
%             Lower = repmat(Global.lower_domain(1,:),N,1);
%             Upper = repmat(Global.lower_domain(2,:),N,1);
%     end
    Site  = rand(N,D) < proM;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
                  
    Offspring       = min(max(Offspring,Lower),Upper);
%     if calObj
%         Offspring = INDIVIDUAL(Offspring,type);
%     end
end
