classdef DS4 < PROBLEM
    methods
        function obj = DS4()
            obj.Global.M = [2 2];
            if isempty(obj.Global.D)
                K = 5;
                L = 4;
                obj.Global.D = [1,K+L];
            else
                if obj.Global.D(1) == 1
                    K = ceil(obj.Global.D(2)/2);
                    L = obj.Global.D(2)-K;
                else
                    assert(obj.Global.D(1) == obj.Global.D(2),'PROBLEM:Error','The upper level dimension must equal to 1 for %s',class(obj));
                end
            end
            obj.Parameter = table(K,L);
            obj.Global.upper_domain = [ones(1,obj.Global.D(1));2*ones(1,obj.Global.D(1))];
            obj.Global.lower_domain = [0,-(K+L)*ones(1,obj.Global.D(2)-1);1,(K+L)*ones(1,obj.Global.D(2)-1)];
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
            L = obj.Parameter.('L');
            
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
                    x1 = x(:,1);
                    temp = sum(x(:,2:K).^2,2)+1;
                    
                    objs(:,1) = (1-x1).*temp.*y1;
                    objs(:,2) = x1.*temp.*y1;
                case 'lower'
                    y1 = y(:,1);
                    x1 = x(:,1);
                    temp = sum(x(:,K+1:K+L).^2,2)+1;
                    objs(:,1) = (1-x1).*temp.*y1;
                    objs(:,2) = x1.*temp.*y1;
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
                    y1 = y(:,1);
                    x1 = x(:,1);
                    cons =-((1-x1).*y1+1/2*x1.*y1-1);
                case 'lower'
                    cons = zeros(N,1);
            end
        end
        
        function P = PF(obj)
			file = fullfile('Problems','DS','DS4.mat');
            try
                load(file);
            catch    
                N = 1e3;
                P(:,1) = (0:1/(N-1):1)';
                P(:,2) = -2*P(:,1)+2;
                save(file,'P');
            end
        end
        
        function P = lower_PF(obj,y)
            N = 1e3;
            y1 = y(:,1);
            t=(0:1/(N-1):y1)';
            P(:,1) = t;
            P(:,2) = y1-t;
        end
    end
    
end