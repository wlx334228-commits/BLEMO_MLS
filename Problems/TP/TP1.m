classdef TP1 < PROBLEM
    methods
        function obj = TP1()
            obj.Global.M = [2;2];
            if isempty(obj.Global.D)
                obj.Global.D = [1;2];
            end
            
            obj.Global.upper_domain = [zeros(1,obj.Global.D(1));ones(1,obj.Global.D(1))];
            obj.Global.lower_domain = [-ones(1,obj.Global.D(2));ones(1,obj.Global.D(2))];
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
        
        function PopObj = CalObj(obj,Population,type)
            
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
                    PopObj(:,1) = x(:,1)-y;
                    PopObj(:,2) = x(:,2);
                case 'lower'
                    PopObj(:,1) = x(:,1);
                    PopObj(:,2) = x(:,2);
            end
        end
        
        function PopCon = CalCon(obj,Population,type)
            
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
                    PopCon = -(1+x(:,1)+x(:,2));
                case 'lower'
                    PopCon = x(:,1).^2+ x(:,2).^2-y.^2;
            end
        end
        
        function P = PF(obj)
			file = fullfile('Problems','TP','TP1.mat');
            try
                load(file);
            catch
                N = 1e6;
                y = (1:1/(N-1):sqrt(2))';
                x2 = [-1/2+1/4*sqrt(4*y.^2-4);-1/2-1/4*sqrt(4*y.^2-4)];
                x1 = -1-x2;
                P(:,1) = x1-1/sqrt(2)*[y;y];
                P(:,2) = x2;
                Z=min(P,[],1);
                W = UniformPoint(1e3,2);
                Next = false(1,size(P,1));
                for i=1:1e3
                    index = find(~Next);
                    [~,pos] = min(max((P(index,:)-Z)./repmat(W(i,:),length(index),1),[],2));
                    Next(index(pos)) = true;
                end
                P = P(Next,:);
                save(file,'P');
            end
        end
        
        function P = lower_PF(obj,y)
            N = 1e6;
            t=(2:1/(N-1):3)';
            P(:,1) = y*cos(pi/2*t);
            P(:,2) = y*sin(pi/2*t);
        end
    end
    
end