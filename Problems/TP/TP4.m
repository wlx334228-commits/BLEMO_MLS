classdef TP4 < PROBLEM
    methods
        function obj = TP4()
            obj.Global.M = [2;2];
            if isempty(obj.Global.D)
                obj.Global.D = [2;3];
            end
            obj.Global.upper_domain = [zeros(1,obj.Global.D(1));1e4*ones(1,obj.Global.D(1))];
            obj.Global.lower_domain = [zeros(1,obj.Global.D(2));1e4*ones(1,obj.Global.D(2))];
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
                    objs(:,1) = -(y(:,1)+9*y(:,2)+ 10*x(:,1)+x(:,2)+3*x(:,3));
                    objs(:,2) = -(9*y(:,1)+2*y(:,2)+ 2*x(:,1)+7*x(:,2)+4*x(:,3));
                case 'lower'
                    objs(:,1) = -(4*y(:,1)+6*y(:,2)+ 7*x(:,1)+4*x(:,2)+8*x(:,3));
                    objs(:,2) = -(6*y(:,1)+4*y(:,2)+ 8*x(:,1)+7*x(:,2)+4*x(:,3));
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
                    cons(:,1) =3*y(:,1)+9*y(:,2)+ 9*x(:,1)+5*x(:,2)+3*x(:,3)-1039;
                    cons(:,2) =-4*y(:,1)-y(:,2)+ 3*x(:,1)-3*x(:,2)+2*x(:,3)-94;
                case 'lower'
                    cons(:,1) = 3*y(:,1)-9*y(:,2)-9*x(:,1)-4*x(:,2)-61;
                    cons(:,2) = 5*y(:,1)+9*y(:,2)+ 10*x(:,1)-x(:,2)-2*x(:,3)-924;
                    cons(:,3) = 3*y(:,1)-3*y(:,2)+x(:,2)+5*x(:,3)-420;
            end
        end
        
        function objs = PF(obj)
            y = [146.2955,28.9394];
            x = [0,67.9318,0];
            objs(:,1) = -(y(:,1)+9*y(:,2)+ 10*x(:,1)+x(:,2)+3*x(:,3));
            objs(:,2) = -(9*y(:,1)+2*y(:,2)+ 2*x(:,1)+7*x(:,2)+4*x(:,3));
%             try
%                 load('Problems\TP\TP2.mat');
%             catch
%                 N = 1e6;
%                 y = (1/2:1/(N-1):1)';
%                 x1=y;
%                 P(:,1) = (x1-1).^2 + y.^2;
%                 P(:,2) = (x1-1).^2 + (y-1).^2;
%                 Z=max(P,[],1);
%                 W = UniformPoint(1e3,2);
%                 Next = false(1,size(P,1));
%                 for i=1:1e3
%                     index = find(~Next);
%                     [~,pos] = min(max((Z-P(index,:))./repmat(W(i,:),length(index),1),[],2));
%                     Next(index(pos)) = true;
%                 end
%                 P = P(Next,:);
%                 save('Problems\TP\TP2.mat','P');
%             end
        end
        
        function P = lower_PF(obj,y)
            P = [];
%             N = 1e6;
%             x1=(0:1/(N-1):y)';
%             P(:,1) = x1.^2;
%             P(:,2) = (x1-y).^2;
        end
        
    end
    
end