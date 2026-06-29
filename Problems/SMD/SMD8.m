classdef SMD8 < PROBLEM 
    methods
        function obj = SMD8()
            eps = 1e-6;
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [2 3];
            end
    
            r = floor(obj.Global.D(1)/2);
            p = obj.Global.D(1) - r;
            q = obj.Global.D(2) - r;
            PS = {zeros(1,obj.Global.D(1)),[ones(1,q),zeros(1,r)]};
            obj.Parameter = table(r,p,q,PS);
            
            obj.Global.upper_domain = [repmat(-5,1,p+r);repmat(10,1,p+r)];
            obj.Global.lower_domain = [repmat(-5,1,q+r);repmat(10,1,q+r)];
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
            p = obj.Parameter.('p');
            r = obj.Parameter.('r');
            q = obj.Parameter.('q');
            
            if isa(Population,'INDIVIDUAL')
                xu = Population.upper_decs;
                xl = Population.lower_decs;
            elseif iscell(Population)
                xu = Population{1};
                xl = Population{2};
            end
            
            xu1 = xu(:,1:p);
            xu2 = xu(:,p+1:p+r);
            xl1 = xl(:,1:q);
            xl2 = xl(:,q+1:q+r);
            
            switch type
                case 'upper'
                    term2 = 0;
                    for i=1:q-1
                        term2 = term2 + (xl1(:,i+1) - xl1(:,i).^2).^2 + (xl1(:,i) - 1).^2;
                    end
                    objs = 20+exp(1)-20*exp(-0.2*sqrt(1/p*sum((xu1).^2,2))) - exp(1/p*sum(cos(2*pi*xu1),2))  ...
                        - term2 ...
                        + sum((xu2).^2,2) - sum((xu2 - xl2.^3).^2,2);
                case 'lower'
                    term2 = 0;
                    for i=1:q-1
                        term2 = term2 + (xl1(:,i+1) - xl1(:,i).^2).^2 + (xl1(:,i) - 1).^2;
                    end
                    objs(:,1) =  sum(abs(xu1),2) ...
                        + term2 ...
                        + sum((xu2 - xl2.^3).^2,2);
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
            PS = obj.Parameter.('PS');
            P = CalObj(obj,PS,'upper');
        end
        
        function P = lower_PF(obj,upper_decs)
            p = obj.Parameter.('p');
            r = obj.Parameter.('r');
            q = obj.Parameter.('q');
            N = size(upper_decs,1);
            lower_decs = [zeros(N,q),upper_decs(:,p+1:p+r).^(1/3)];
            P = CalObj(Obj,{upper_decs,lower_decs},'lower');
        end
    
    end
end