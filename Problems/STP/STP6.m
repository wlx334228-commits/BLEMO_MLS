classdef STP6 < PROBLEM 
    methods
        function obj = STP6()
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [1 2];
            end
            
            obj.Global.upper_domain = [0;2];
            obj.Global.lower_domain = [0 0;2 2];
            
            PS = {[1.888],[0.888,0]};
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
            y1 = xl(:,1);
            y2 = xl(:,2);
            
            switch type
                case 'upper'
                    objs = (x1-1).^2 + 2*y1 - 2*x1;
                case 'lower'
                    objs(:,1) = (2*y1-4).^2 + (2*y2-1).^2 + x1*y1;
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
            y1 = xl(:,1);
            y2 = xl(:,2);
            
            switch type
                case 'upper'
                    cons(:,1) =zeros(size(x1,1),1);
                case 'lower'
                    cons(:,1) =  4*x1+5*y1+4*y2-12;
                    cons(:,2) =  4*y2-4*x1-5*y1+4;
                    cons(:,3) =  4*x1-4*y1+5*y2-4;
                    cons(:,3) =  4*y1-4*x1+5*y2-4;
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