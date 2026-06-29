function Offspring = PCX(indexParent,Parent1,Parent2,type,Parameter)
% PCX

    %% Parameter setting
    if nargin > 3
        [CR,proM,disM] = deal(Parameter{:});
    else
        [CR,proM,disM] = deal(0.5,1,20);
    end
    if isa(indexParent,'INDIVIDUAL') && isa(Parent1,'INDIVIDUAL') && isa(Parent2,'INDIVIDUAL')
        calObj = true;
        switch type
            case 'upper'
                indexParent = indexParent.upper_dec;
                Parent1 = Parent1.upper_decs;
                Parent2 = Parent2.upper_decs;
                g = mean([indexParent;Parent1;Parent2]);
            case 'lower'
                indexParent = indexParent.lower_dec;
                Parent1 = Parent1.lower_decs;
                Parent2 = Parent2.lower_decs;
                g = mean([indexParent;Parent1;Parent2]);
        end
    else
        calObj = false;
    end
    
    [N,D]  = size(Parent1);
    epsilon = 1e-10;
    
    sigma1 = 0.1;
    if mean(abs(indexParent - g)) <= epsilon
        sigma2 = 0;
    else
        sigma2 = 1/mean(abs(indexParent - g));
    end
    
    Site = rand(N,D) < CR;
    Offspring = repmat(indexParent,N,1);
    g = repmat(g,N,1);
    Offspring(Site) = Offspring(Site) + randn*sigma1*(Offspring(Site)-g(Site)) + randn*sigma2*(Parent2(Site) - Parent1(Site))/2;
    
    %% Polynomial mutation
    
    Global=GLOBAL.GetObj();
    switch type
        case 'upper'
            Lower = repmat(Global.upper_domain(1,:),N,1);
            Upper = repmat(Global.upper_domain(2,:),N,1);
        case 'lower'
            Lower = repmat(Global.lower_domain(1,:),N,1);
            Upper = repmat(Global.lower_domain(2,:),N,1);
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
