classdef DS5 < PROBLEM
    methods
        function obj = DS5()
            alpha = 5;
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
            obj.Parameter = table(K,L,alpha);
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
            alpha = obj.Parameter.('alpha');
            
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
                    cons =-((1-x1).*y1+x1.*y1-2+1/alpha*ceil(alpha*(1-x1).*y1));
                case 'lower'
                    cons = zeros(N,1);
            end
        end
        
        function P = PF(obj)
			file = fullfile('Problems','DS','DS5.mat');
            try
                load(file);
            catch    
                N = 2e3;
                y1=[1,1.2,1.4,1.6,1.8];
                lower = cell2mat(cellfun(@(x) 2*(1-1/x), num2cell(y1), 'UniformOutput',false));
                upper = cell2mat(cellfun(@(x) 2*(1-0.9/x), num2cell(y1), 'UniformOutput',false));
                P = [];
                for i=1:length(y1)
                    t = (lower(i):1/(N-1):upper(i))';
                    P = cat(1,P,[(1-t).*y1(i),t.*y1(i)]);
                end
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
            y1 = y(:,1);
            t=(0:1/(N-1):y1)';
            P(:,1) = t;
            P(:,2) = y1-t;
        end
    end
    
end