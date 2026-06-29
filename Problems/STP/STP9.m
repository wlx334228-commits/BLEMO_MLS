classdef STP9 < PROBLEM 
    methods
        function obj = STP9()
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [5 5];
            end
            
            obj.Global.upper_domain = [-ones(1,obj.Global.D(1));ones(1,obj.Global.D(1))];
            obj.Global.lower_domain = [-pi*ones(1,obj.Global.D(2));pi*ones(1,obj.Global.D(2))];
            
        end
        
        function Decs = Decs(obj,Decs,type)
            switch type
                case 'upper'
                    domain = obj.Global.upper_domain;
                case 'lower'
                    domain = obj.Global.lower_domain;
                case 'bilevel'
                    domain = [obj.Global.upper_domain,obj.Global.lower_domain];
            end
            Lower = repmat(domain(1,:),length(obj),1);
            Upper = repmat(domain(2,:),length(obj),1);
            Decs  = max(min(Decs,Upper),Lower);
        end 
        
        function objs = CalObj(obj,Population,type)
            
            if isa(Population,'INDIVIDUAL')
                x = Population.upper_decs;
                y = Population.lower_decs;
            elseif iscell(Population)
                x = Population{1};
                y = Population{2};
            end
            
            switch type
                case 'upper'
                    objs = sum(abs(x-1),2)+sum(abs(y),2);
                case 'lower'
                    exponent = 1+1/4000*sum(y.^2,2)-prod(cos(y./sqrt(1:obj.Global.D(2))),2);
                    exponent = exponent*sum(x.^2,2);
                    objs(:,1) = exp(exponent);
            end
        end
        
        function cons = CalCon(obj,Population,type)
            
            if isa(Population,'INDIVIDUAL')
                xu = Population.upper_decs;
                xl = Population.lower_decs;
            elseif iscell(Population)
                xu = Population{1};
                xl = Population{2};
            end
            
            switch type
                case 'upper'
                    cons =zeros(size(xu,1),1);
                    
                case 'lower'
                    cons(:,1) =  zeros(size(xl,1),1);
            end
        end
        
        function P = PF(obj,type)
            P = 0;
            switch type
                case 'bilevel'
                    P=[0,1];
            end
        end
        
        function P = lower_PF(obj,upper_decs)
            P= [];
        end
    
    end
end