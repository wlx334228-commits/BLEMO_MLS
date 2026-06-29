classdef DS2_2 < PROBLEM
    methods
        function obj = DS2_2()
            r = 0.25;
            gamma = 4;
            tao = -1;
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
            obj.Parameter = table(K,r,gamma,tao);
            obj.Global.upper_domain = [0.001,-K*ones(1,obj.Global.D(1)-1);K,K*ones(1,obj.Global.D(1)-1)];
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
                    v1=zeros(N,1);
                    v2=zeros(N,1);
                    y1 = y(:,1);
                    v1(y1>1) = y1(y1>1)-(1-cos(0.2*pi));
                    v1(~(y1>1)) = cos(0.2*pi)*y1(~(y1>1))+sin(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1(~(y1>1)))));
                    v2(y1>1) = 0.1*(y1(y1>1)-1)-sin(0.2*pi);
                    v2(~(y1>1)) = -sin(0.2*pi)*y1(~(y1>1))+cos(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1(~(y1>1)))));
                    
                    objs(:,1) = v1+sum(y(:,2:end).^2+10*(1-cos(pi/K*y(:,2:end))),2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*cos(gamma*pi/2*x(:,1)./y(:,1));
                    objs(:,2) = v2+sum(y(:,2:end).^2+10*(1-cos(pi/K*y(:,2:end))),2)+tao*sum((x(:,2:end)-...
                                 y(:,2:end)).^2,2)-r*sin(gamma*pi/2*x(:,1)./y(:,1));
                case 'lower'
                    objs(:,1) = x(:,1).^2 + sum((x(:,2:end)-y(:,2:end)).^2,2);
                    objs(:,2) = sum((x-y).^2.*repmat(1:K,N,1),2);
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
		
			file = fullfile('Problems','DS','DS2.mat');
            try
                load(file);
            catch
                N = 1e3;
                r = obj.Parameter.('r');
                y1=[0.001 0.2 0.4 0.6 0.8 1];
                v1 = cos(0.2*pi)*y1+sin(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1)));
                v2 = -sin(0.2*pi)*y1+cos(0.2*pi)*sqrt(abs(0.02*sin(5*pi*y1)));
                t=(2:1/(N-1):3)';
                P(:,1) = reshape(cell2mat(cellfun(@(x) x+r*cos(pi/2*t), num2cell(v1), 'UniformOutput',false)),[],1);
                P(:,2) = reshape(cell2mat(cellfun(@(x) x+r*sin(pi/2*t), num2cell(v2), 'UniformOutput',false)),[],1);
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
%             r = obj.Parameter.('r');
            t=(0:1/(N-1):y(:,1))';
            P(:,1) = t.^2;
            P(:,2) = (t-y(:,1)).^2;
        end
    end
    
end