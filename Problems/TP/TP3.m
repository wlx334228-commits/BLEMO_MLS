classdef TP3 < PROBLEM
    methods
        function obj = TP3()
            obj.Global.M = [2;2];
            if isempty(obj.Global.D)
                obj.Global.D = [1;2];
            end
            
            obj.Global.upper_domain = [zeros(1,obj.Global.D(1));10*ones(1,obj.Global.D(1))];
            obj.Global.lower_domain = [zeros(1,obj.Global.D(2));10*ones(1,obj.Global.D(2))];
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
                N = length(Population);
                y = Population.upper_decs;
                x = Population.lower_decs;
            elseif iscell(Population)
                y = Population{1};
                x = Population{2};
                N = size(x,1);
            end
            
            switch type
                case 'upper'
                    objs(:,1) = x(:,1)+x(:,2).^2 + y(:,1) + sin(x(:,1)+y(:,1)).^2;
                    objs(:,2) = cos(x(:,2)).*(0.1+y(:,1)).*exp(-x(:,1)./(0.1+x(:,2)));
                case 'lower'
                    objs(:,1) = ((x(:,1)-2).^2+(x(:,2)-2).^2)/4 + (x(:,2).*y(:,1)+(5-y(:,1)).^2)/16+sin(x(:,2)/10);
                    objs(:,2) = (x(:,1).^2+(x(:,2)-6).^4-2*x(:,1).*y(:,1)-(5-y(:,1)).^2)/80;
            end
        end
        
        function cons = CalCon(obj,Population,type)
            if isa(Population,'INDIVIDUAL')
                N = length(Population);
                y = Population.upper_decs;
                x = Population.lower_decs;
            elseif iscell(Population)
                y = Population{1};
                x = Population{2};
                N = size(x,1);
            end
            switch type
                case 'upper'
                    cons =(x(:,1)-0.5).^2+(x(:,2)-5).^2+(y(:,1)-5).^2-16;
                case 'lower'
                    cons(:,1) = -x(:,2)+x(:,1).^2;
                    cons(:,2) = 5*x(:,1).^2+x(:,2)-10;
                    cons(:,3) = x(:,2)+y(:,1)/6-5;
                    cons(:,4) = -x(:,1);
            end
        end
        
        function P = PF(obj)
            P=[];
        end
        
        function P = lower_PF(obj,y)
            P = [];
        end
        
    end
    
end