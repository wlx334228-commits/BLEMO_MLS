classdef STP3 < PROBLEM 
    methods
        function obj = STP3()
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [2 2];
            end
            
            obj.Global.upper_domain = [0 0;10 10];
            obj.Global.lower_domain = [0 0;10 10];
            
            PS = {[0,2],[1.8750,0.9062]};
            obj.Parameter = table(PS);
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
                    objs = -x1.^2-3*x2.^2-4*y1+y2.^2;
                case 'lower'
                    objs(:,1) = 2*x1.^2+y1.^2-5*y2;
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
                    cons(:,1) =x1.^2+2*x2-4;
                case 'lower'
                    cons(:,1) =  -3-x1.^2+2*x1-x2.^2+2*y1-y2;
                    cons(:,2) =  4-x2-3*y1+4*y2;
            end
        end
        
        function P = PF(obj,type)
            PS = obj.Parameter.('PS');
            P = CalObj(obj,PS,'upper');
            if nargin>1
                switch type
                    case 'bilevel'
                        P=[P,CalObj(obj,PS,'lower')];
                end
            end
        end
        
        function P = lower_PF(obj,upper_decs)
            P= [];
        end
    
    end
end