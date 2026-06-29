classdef DS1_1 < PROBLEM
    methods
        function obj = DS1_1()   
            r = 0.1;     
            alpha = 1;   
            gamma = 1;
            tao = 1;           
%             tao = -1;
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
            obj.Parameter = table(K,r,alpha,gamma,tao);
            obj.Global.upper_domain = [1,-K*ones(1,obj.Global.D(1)-1);4,K*ones(1,obj.Global.D(1)-1)];
            obj.Global.lower_domain = [-K*ones(1,obj.Global.D(2));K*ones(1,obj.Global.D(2))];
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
            K = obj.Parameter.('K');
            r = obj.Parameter.('r');
            alpha = obj.Parameter.('alpha');
            gamma = obj.Parameter.('gamma');
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
                    temp = repmat(cell2mat(cellfun(@(x) (x-1)/2, {2:K}, 'UniformOutput',false)),N,1);
                    objs(:,1) = (1+r-cos(alpha*pi*y(:,1)))+sum((y(:,2:end)-temp).^2,2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*cos(gamma*pi/2*x(:,1)./y(:,1));
                    objs(:,2) = (1+r-sin(alpha*pi*y(:,1)))+sum((y(:,2:end)-temp).^2,2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*sin(gamma*pi/2*x(:,1)./y(:,1));
                case 'lower'
                    objs(:,1) = x(:,1).^2 + sum((x(:,2:end)-y(:,2:end)).^2,2)+10*sum(1-cos(pi/K*(x(:,2:end)-y(:,2:end))),2);
                    objs(:,2) = sum((x-y).^2,2)+10*sum(abs(sin(pi/K*(x(:,2:end)-y(:,2:end)))),2);
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
                    cons =zeros(N,1);
                case 'lower'
                    cons = zeros(N,1);
            end
        end
        
        function P = PF(obj)
			file = fullfile('Problems','DS','DS1.mat');
            try
                load(file);
            catch
                N = 1e3;
                r = obj.Parameter.('r');
                t=(2:1/(N-1):3)';
                P(:,1) = 1+r+(1+r)*cos(pi/2*t);
                P(:,2) = 1+r+(1+r)*sin(pi/2*t);
                save(file,'P');
            end
        end
        
        function P = lower_PF(obj,y)
            N = 1e6;
%             r = obj.Parameter.('r');
            t=(0:1/(N-1):y(:,1))';
            P(:,1) = t.^2;
            P(:,2) = (t-y(:,1)).^2;
        end
    end
    
end