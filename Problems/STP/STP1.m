classdef STP1 < PROBLEM 
    methods
        function obj = STP1()
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [2 2];
            end
            
            obj.Global.upper_domain = [-30 -30;30 15];
            obj.Global.lower_domain = [0 0;10 10];
            
            PS = {[20,5],[10,5]};
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
            
            xu1 = xu(:,1);
            xu2 = xu(:,2);
            xl1 = xl(:,1);
            xl2 = xl(:,2);
            
            switch type
                case 'upper'
                    objs = (xu1-30).^2 + (xu2-20).^2 - 20*xl1 + 20*xl2;
                case 'lower'
                    objs(:,1) = (xu1 - xl1).^2 + (xu2 - xl2).^2;
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
            
            xu1 = xu(:,1);
            xu2 = xu(:,2);
            xl1 = xl(:,1);
            xl2 = xl(:,2);
            
            switch type
                case 'upper'
                    cons(:,1) =30 - xu1 - 2*xu2;
                    cons(:,2) =xu1 + xu2 - 25;
                case 'lower'
                    cons = zeros(size(xu,1),1);
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