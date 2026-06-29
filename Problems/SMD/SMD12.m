classdef SMD12 < PROBLEM 
    methods
        function obj = SMD12()
            eps = 1e-6;
            obj.Global.M = [1 1];
            if isempty(obj.Global.D)
                obj.Global.D = [2 3];
            end
    
            r = floor(obj.Global.D(1)/2);
            p = obj.Global.D(1) - r;
            q = obj.Global.D(2) - r;
            PS = {repmat(1/sqrt(p+r-1),1,obj.Global.D(1)),[repmat(1/sqrt(q-1),1,q),repmat(atan(1/sqrt(p+r-1)-1/sqrt(r)),1,r)]};
            obj.Parameter = table(r,p,q,PS);
            
            obj.Global.upper_domain = [repmat(-5,1,p),-14.10*ones(1,r);repmat(10,1,p),14.10*ones(1,r)];
            obj.Global.lower_domain = [repmat(-5,1,q),repmat(-1.5+eps,1,r);repmat(10,1,q),repmat(1.5-eps,1,r)];
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
                    a = 2*ones(size(xu1));
                    c = 2*ones(size(xu2));
                    objs = sum((xu1 - a).^2,2) ...
                        + sum((xl1).^2,2) ...
                        + sum((xu2 - c).^2,2) + sum(tan(abs(xl2)),2) - sum((xu2 - tan(xl2)).^2,2);
                    
                case 'lower'
                    b = 2*ones(size(xl1));
                    objs(:,1) = sum((xu1).^2,2) ...
                        + sum((xl1 - b).^2,2) ...
                        + sum((xu2 - tan(xl2)).^2,2);
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
                        cons(:,i)= -xu2(:,i) + tan(xl2(:,i));
                    end
                    
                    for i=1:p
                        cons(:,r+i) = -(xu1(:,i) + xu1(:,i).^3 - sum(xu1.^3,2) - sum(xu2.^3,2));
                    end
                    
                    for i=1:r
                        cons(:,r+p+i) = -(xu2(:,i) + xu2(:,i).^3 - sum(xu2.^3,2) - sum(xu1.^3,2));
                    end
                case 'lower'
                    for i=1:q
                        cons(:,i)=-(xl1(:,i) + xl1(:,i).^3 - sum(xl1.^3,2));
                    end
                    
                    cons(:,q+1)=-sum((xu2 - tan(xl2)).^2,2) + 1;
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