classdef SMD11 < PROBLEM 
    methods
        function obj = SMD11()
            eps = 1e-6;
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [2 3];
            end
    
            r = floor(obj.Global.D(1)/2);
            p = obj.Global.D(1) - r;
            q = obj.Global.D(2) - r;
            PS = {zeros(1,obj.Global.D(1)),[zeros(1,q),repmat(exp(-1/sqrt(r)),1,r)]};
            obj.Parameter = table(r,p,q,PS);
            
            obj.Global.upper_domain = [repmat(-5,1,p),-ones(1,r);repmat(10,1,p),ones(1,r)];
            obj.Global.lower_domain = [repmat(-5,1,q),repmat(1/exp(1),1,r);repmat(10,1,q),repmat(exp(1),1,r)];
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
                    objs = sum((xu1).^2,2) ...
                               - sum((xl1).^2,2) ...
                               + sum((xu2).^2,2) - sum((xu2 - log(xl2)).^2,2);
                case 'lower'
                    objs(:,1) = sum((xu1).^2,2) ...
                               + sum((xl1).^2,2) ...
                               + sum((xu2 - log(xl2)).^2,2);
            end
        end
        
        function cons = CalCon(obj,Population,type)
            
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
                    for i=1:r
                        cons(:,i) = -(xu2(:,i) - 1/sqrt(r) - log(xl2(:,i)));
                    end
                case 'lower'
                    cons = -(sum((xu2 - log(xl2)).^2,2) - 1);
            end
        end
        
        function P = PF(obj)
            PS = obj.Parameter.('PS');
            P = CalObj(obj,PS,'upper');
        end
        
        function P = lower_PF(obj,upper_decs)
            P = NaN;
        end
    
    end
end