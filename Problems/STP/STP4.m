classdef STP4 < PROBLEM 
    methods
        function obj = STP4()
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [2 3];
            end
            
            obj.Global.upper_domain = [0 0;1 1];
            obj.Global.lower_domain = [0 0 0;1 1 1];
            
%             PS = {[0.29,0.70],[0,0.27,0.27]};
%             obj.Parameter = table(PS);
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
            y3 = xl(:,3);
            
            switch type
                case 'upper'
                    objs = -8*x1-4*x2+4*y1-40*y2-4*y3;
                case 'lower'
                    objs(:,1) = x1+2*x2+y1+y2+2*y3;
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
            y3 = xl(:,3);
            
            switch type
                case 'upper'
                    cons(:,1) =zeros(size(x1,1),1);
                case 'lower'
                    cons(:,1) =  y2+y3-y1-1;
                    cons(:,2) =  2*x1-y1+2*y2-0.5*y3-1;
                    cons(:,3) =  2*x2+2*y1-y2-0.5*y3-1;
            end
        end
        
        function P = PF(obj,type)
%             PS = obj.Parameter.('PS');
            P = -29.2;
            if nargin>1
                switch type
                    case 'bilevel'
                        P=[-29.2,3.2];
                end
            end
        end
        
        function P = lower_PF(obj,upper_decs)
            P= [];
        end
    
    end
end