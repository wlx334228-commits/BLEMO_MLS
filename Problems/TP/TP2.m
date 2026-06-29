classdef TP2 < PROBLEM
    methods
        function obj = TP2()
            obj.Global.M = [2;2];
            if isempty(obj.Global.D)
                obj.Global.D = [1;14];
            end
            
            obj.Global.upper_domain = [-ones(1,obj.Global.D(1));2*ones(1,obj.Global.D(1))];
            obj.Global.lower_domain = [-ones(1,obj.Global.D(2));2*ones(1,obj.Global.D(2))];
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
                    objs(:,1) = (x(:,1)-1).^2 + sum(x(:,2:end).^2,2)+y.^2;
                    objs(:,2) = (x(:,1)-1).^2 + sum(x(:,2:end).^2,2)+(y-1).^2;
                case 'lower'
                    objs(:,1) = x(:,1).^2 + sum(x(:,2:end).^2,2);
                    objs(:,2) = (x(:,1)-y).^2 + sum(x(:,2:end).^2,2);
            end
        end
        
        function cons = CalCon(obj,Population,type)
            N = length(Population);
            
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
                    cons =zeros(N,1);
                case 'lower'
                    cons = zeros(N,1);
            end
        end
        
        function P = PF(obj)
			file = fullfile('Problems','TP','TP2.mat');
            try
                load(file);
            catch
                N = 1e6;
                y = (1/2:1/(N-1):1)';
                x1=y;
                P(:,1) = (x1-1).^2 + y.^2;
                P(:,2) = (x1-1).^2 + (y-1).^2;
                Z=max(P,[],1);
                W = UniformPoint(1e3,2);
                Next = false(1,size(P,1));
                for i=1:1e3
                    index = find(~Next);
                    [~,pos] = min(max((Z-P(index,:))./repmat(W(i,:),length(index),1),[],2));
                    Next(index(pos)) = true;
                end
                P = P(Next,:);
                save(file,'P');
            end
        end
        
        function P = lower_PF(obj,y)
            N = 1e6;
            x1=(0:1/(N-1):y)';
            P(:,1) = x1.^2;
            P(:,2) = (x1-y).^2;
        end
        
    end
    
end