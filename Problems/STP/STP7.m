classdef STP7 < PROBLEM 
    methods
        function obj = STP7()
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [2 2];
            end
            
            obj.Global.upper_domain = [0 0;10 10];
            obj.Global.lower_domain = [0 0;1 10];
            
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
                xu = Population.upper_decs;
                xl = Population.lower_decs;
            elseif iscell(Population)
                xu = Population{1};
                xl = Population{2};
            end
            
            x1 = xu(:,1);
            x2 = xu(:,2);
            y1 = xl(:,1);
            y2 = xl(:,2);
            
            switch type
                case 'upper'
                    objs = -(x1+y1).*(x2+y2)./(1+x1.*y1+x2.*y2);
                case 'lower'
                    objs(:,1) = (x1+y1).*(x2+y2)./(1+x1.*y1+x2.*y2);
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
            
            x1 = xu(:,1);
            x2 = xu(:,2);
            y1 = xl(:,1);
            y2 = xl(:,2);
            
            switch type
                case 'upper'
                    cons =[x1.^2+x2.^2 - 100,x1-x2];
                    
                case 'lower'
                    cons(:,1) =  y1-x1;
                    cons(:,2) =  y2-x2;
            end
        end
        
        function P = PF(obj,type)
            P = -1.96;
            if nargin>1
                switch type
                    case 'bilevel'
                        P=[-1.96,1.96];
                end
            end
        end
        
        function P = lower_PF(obj,upper_decs)
            P= [];
        end
    
    end
end