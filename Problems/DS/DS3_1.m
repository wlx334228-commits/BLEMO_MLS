classdef DS3_1 < PROBLEM
    methods
        function obj = DS3_1()
            r = 0.2;
            tao = 1;
            obj.Global.M = [2 2];
            if isempty(obj.Global.D)
                K = 10;
                obj.Global.D = [K,K];
            else
                if obj.Global.D(1) == obj.Global.D(2)
                    K = obj.Global.D(1);
                else
                    assert(obj.Global.D(1) == obj.Global.D(2),'PROBLEM:Error','The upper level dimension must equal to the lower level dimension for %s',class(obj));
                end
            end
            obj.Parameter = table(K,r,tao);
            obj.Global.upper_domain = [zeros(1,obj.Global.D(1));K*ones(1,obj.Global.D(1))];
            obj.Global.lower_domain = [-K*ones(1,obj.Global.D(2));K*ones(1,obj.Global.D(2))];
        end
        
        function Decs = Decs(obj,Decs,type)
            switch type
                case 'upper'
                    domain = obj.Global.upper_domain;
                    Lower = repmat(domain(1,:),length(obj),1);
                    Upper = repmat(domain(2,:),length(obj),1);
                    Decs  = max(min(Decs,Upper),Lower);
                    Decs(:,1) = round(Decs(:,1),1);
                case 'lower'
                    domain = obj.Global.lower_domain;
                    Lower = repmat(domain(1,:),length(obj),1);
                    Upper = repmat(domain(2,:),length(obj),1);
                    Decs  = max(min(Decs,Upper),Lower);
                case 'bilevel'
                    domain = [obj.Global.upper_domain,obj.Global.lower_domain];
                    Lower = repmat(domain(1,:),length(obj),1);
                    Upper = repmat(domain(2,:),length(obj),1);
                    Decs  = max(min(Decs,Upper),Lower);
                    Decs(:,1) = round(Decs(:,1),1);
            end
            
        end
        
        function objs = CalObj(obj,Population,type)
            K = obj.Parameter.('K');
            tao = obj.Parameter.('tao');
            
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
                    y1 = y(:,1);
                    y2 = y(:,2);
                    R = 0.1+0.15*abs(sin(2*pi*(y1-0.1)));
                    temp = repmat(cell2mat(cellfun(@(x) x/2, {3:K}, 'UniformOutput',false)),N,1);
                    objs(:,1) = y1+sum((y(:,3:end)-temp).^2,2)+tao*sum((x(:,3:end)-...
                                 y(:,3:end)).^2,2)-R.*cos(4*atan((y(:,2)-x(:,2))./(y(:,1)-x(:,1)+eps)));
                    objs(:,2) = y2+sum((y(:,3:end)-temp).^2,2)+tao*sum((x(:,3:end)-...
                                 y(:,3:end)).^2,2)-R.*sin(4*atan((y(:,2)-x(:,2))./(y(:,1)-x(:,1)+eps)));
                case 'lower'
                    objs(:,1) = x(:,1) + sum((x(:,3:end)-y(:,3:end)).^2,2);
                    objs(:,2) = x(:,2) + sum((x(:,3:end)-y(:,3:end)).^2,2);
            end
        end
        
        function cons = CalCon(obj,Population,type)
            r = obj.Parameter.('r');
            
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
                    cons =(1-y(:,1).^2)-y(:,2);
                case 'lower'
                    cons = (y(:,1)-x(:,1)).^2+(y(:,2)-x(:,2)).^2-r^2;
            end
        end
        
        function P = PF(obj)
			file = fullfile('Problems','DS','DS3.mat');
            try
                load(file);
            catch
                N = 1000;
                y1 = (0:0.1:1.3)';
                y2 = max(1-y1.^2,0);
                R = 0.1+0.15*abs(sin(2*pi*(y1-0.1)));
                t=1+0.5*(0:1/(N-1):1)';
                P(:,1) = cell2mat(cellfun(@(x,y) y+x.*cos(pi*t), num2cell(R), num2cell(y1),'UniformOutput',false));
                P(:,2) = cell2mat(cellfun(@(x,y) y+x.*sin(pi*t), num2cell(R), num2cell(y2),'UniformOutput',false));
                Front = NDSort(P,1);
                P = P(Front==1,:);
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
            r = obj.Parameter.('r');
            t=(2:1/(N-1):3)';
            P(:,1) = y(:,1) + r*cos(pi/2*t);
            P(:,2) = y(:,2) + r*sin(pi/2*t);
        end
    end
    
end