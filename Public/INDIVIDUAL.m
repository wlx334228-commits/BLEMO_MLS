classdef INDIVIDUAL < handle
    properties(SetAccess = private)
        upper_dec;        % Decision variables of the individual
        lower_dec
        upper_obj;        % Objective values of the individual
        lower_obj
        upper_con;        % Constraint violations of the individual
        lower_con;
        upper_feasible = false;
        lower_feasible = false;
        upper_CV;
        lower_CV;
        upper_CalObj = false;
        lower_CalObj = false;
        upper_CalCon = false;
        lower_CalCon = false;
        upper_add;
        lower_add;
    end
    
    properties
        label = false;
    end
    methods
        function obj = INDIVIDUAL(Decs,type)
            if nargin > 0
                Global = GLOBAL.GetObj();
                if ~iscell(Decs)
                    obj(1,size(Decs,1)) = INDIVIDUAL;
                    Decs = Global.problem.Decs(Decs,type);
                    switch type
                        case 'upper'
                            for i = 1 : length(obj)
                                obj(i).upper_dec = Decs(i,:);
                            end
                        case 'lower'
                            for i = 1 : length(obj)
                                obj(i).lower_dec = Decs(i,:);
                            end
                    end
                else
                    upper_decs = Global.problem.Decs(Decs{1},'upper');
                    lower_decs = Global.problem.Decs(Decs{2},'lower');
                    if size(upper_decs,1) == size(lower_decs,1)
                        obj(1,size(upper_decs,1)) = INDIVIDUAL;
                        for i = 1 : length(obj)
                            obj(i).upper_dec = upper_decs(i,:);
                            obj(i).lower_dec = lower_decs(i,:);
                        end
                    elseif size(upper_decs,1)==1
                        obj(1,size(lower_decs,1)) = INDIVIDUAL;
                        for i = 1 : length(obj)
                            obj(i).upper_dec = upper_decs;
                            obj(i).lower_dec = lower_decs(i,:);
                        end
                    end
                end
            end
        end
        
        function obj = get_lower_decs(obj,lower_decs)
            obj = INDIVIDUAL({obj.upper_decs, lower_decs});
        end
        
        function obj = update_value(obj,value,type1,type2)
            L = length(obj);
            switch type1
                case 'upper'
                    switch type2
                        case 'obj'
                            for i = 1:L
                                obj(i).upper_obj = value(i,:);
                                obj(i).upper_CalObj = true;
                            end
                        case 'con'
                            for i = 1:L
                                obj(i).upper_con = value(i,:);
                                obj(i).upper_CalCon = true;
                                obj(i).upper_CV = sum(max(0,value(i,:)),2);
                                if obj(i).upper_CV==0
                                    obj(i).upper_feasible = true;
                                end
                            end
                        case 'add'
                            for i = 1:L
                                obj(i).upper_add = value(i,:);
                            end
                    end
                case 'lower'
                    switch type2
                        case 'obj'
                            for i = 1:L
                                obj(i).lower_obj = value(i,:);
                                obj(i).lower_CalObj = true;
                            end
                        case 'con'
                            for i = 1:L
                                obj(i).lower_con = value(i,:);
                                obj(i).lower_CalCon = true;
                                obj(i).lower_CV = sum(max(0,value(i,:)),2);
                                if obj(i).lower_CV==0
                                    obj(i).lower_feasible = true;
                                end
                            end
                        case 'add'
                            for i = 1:L
                                obj(i).lower_add = value(i,:);
                            end
                    end
            end
            
        end
        %% Get the matrix of decision variables of the population
        function value = upper_decs(obj,dim)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            if ~isempty(obj)
                value = cat(1,obj.upper_dec);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = lower_decs(obj,dim)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            if ~isempty(obj)
                value = cat(1,obj.lower_dec);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = decs(obj,dim)
            if ~isempty(obj)
                value = [obj.upper_decs,obj.lower_decs];
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        %% Get the matrix of objective values of the population
        function value = upper_objs(obj,dim)
        %objs - Get the matrix of objective values of the population.
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            if ~isempty(obj)
                value = cat(1,obj.upper_obj);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = lower_objs(obj,dim)
        %objs - Get the matrix of objective values of the population.
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            if ~isempty(obj)
                value = cat(1,obj.lower_obj);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = objs(obj,dim)
            if ~isempty(obj)
                value = [obj.upper_objs,obj.lower_objs];
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        %% Get the matrix of constraint violations of the population
        function value = upper_cons(obj,dim)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            if ~isempty(obj)
                value = cat(1,obj.upper_con);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = lower_cons(obj,dim)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            if ~isempty(obj)
                value = cat(1,obj.lower_con);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = cons(obj,dim)
            if ~isempty(obj)
                value = [obj.upper_cons,obj.lower_cons];
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        %%
        function [value,pop] = bilevel_feasibles(obj,type)
            if ~isempty(obj)
                if nargin > 1
                    switch type
                        case 'true'
                            value = obj.upper_feasibles & obj.lower_feasibles;
                            pop = obj(value);
                        case 'false'
                            value = obj.upper_feasibles & obj.lower_feasibles;
                            pop = obj(~value);
                    end
                else
                    value = obj.upper_feasibles & obj.lower_feasibles;
                    pop = obj(value);
                end
                
            else
                value = [];
                pop = [];
            end
        end
        
        %%
        function [value,pop] = upper_feasibles(obj,type)
            if ~isempty(obj)
                if nargin > 1
                    switch type
                        case 'true'
                            value = cat(2,obj.upper_feasible);
                            pop = obj(value);
                        case 'false'
                            value = cat(2,obj.upper_feasible);
                            pop = obj(~value);
                    end
                else
                    value = cat(2,obj.upper_feasible);
                    pop = obj(value);
                end
                
            else
                value = [];
                pop = [];
            end
        end
        
         %%
        function [value,pop] = lower_feasibles(obj,type)
            if ~isempty(obj)
                if nargin > 1
                    switch type
                        case 'true'
                            value = cat(2,obj.lower_feasible);
                            pop = obj(value);
                        case 'false'
                            value = cat(2,obj.lower_feasible);
                            pop = obj(~value);
                    end
                else
                    value = cat(2,obj.lower_feasible);
                    pop = obj(value);
                end
            else
                value = [];
                pop = [];
            end
        end
        
         %%
        function value = upper_CVs(obj,dim)
            if ~isempty(obj)
                value = cat(1,obj.upper_CV);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
         %%
        function value = lower_CVs(obj,dim)
            if ~isempty(obj)
                value = cat(1,obj.lower_CV);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        %% Get the matrix of additional properties of the population
        function [labels,pop] = labels(obj,type)
            
            if ~isempty(obj)
                if nargin > 1
                    switch type
                        case true
                            labels = cat(2,obj.label);
                            pop = obj(labels);
                        case false
                            labels = cat(2,obj.label);
                            labels = ~labels;
                            pop = obj(labels);
                    end
                else
                    labels = cat(2,obj.label);
                    pop = obj(labels);
                end
                
            else
                pop = [];
                labels = [];
            end
        end
        
        function value = upper_CalObjs(obj)
            
            if ~isempty(obj)
                value = cat(2,obj.upper_CalObj);
            else
                value = [];
            end
        end
        
        function value = upper_CalCons(obj)
            
            if ~isempty(obj)
                value = cat(2,obj.upper_CalCon);
            else
                value = [];
            end
        end
        
        function value = lower_CalObjs(obj,dim)
            
            if ~isempty(obj)
                value = cat(2,obj.lower_CalObj);
                if nargin > 2
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = lower_CalCons(obj,dim)
            
            if ~isempty(obj)
                value = cat(2,obj.lower_CalCon);
                if nargin > 2
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function [value,Best_ind] = bilevel_best(obj)
            index = find(obj.bilevel_feasibles);
            if isempty(index)
                [~,Best_ind] = min(obj.upper_CVs);
                value = obj(Best_ind);
            elseif length(obj(index(1)).upper_obj) > 1
                Best = NDSort(obj(index).upper_objs,1);
                Best_ind = index(Best == 1);
                value = obj(Best_ind);
            else
                [~,Best] = min(obj(index).upper_objs);
                Best_ind = index(Best);
                value = obj(Best_ind);
            end
            
        end
        
        function [value,Best_ind] = upper_best(obj)
            index = find(obj.upper_feasibles);
            if isempty(index)
                [~,Best_ind] = min(obj.upper_CVs);
                value = obj(Best_ind);
            elseif length(obj(index(1)).upper_obj) > 1
                Best = NDSort(obj(index).upper_objs,1);
                Best_ind = index(Best == 1);
                value = obj(Best_ind);
            else
                [~,Best] = min(obj(index).upper_objs);
                Best_ind = index(Best);
                value = obj(Best_ind);
            end
            
        end
        
        function [value,Best_ind] = lower_best(obj)
            index = find(obj.lower_feasibles);
            if isempty(index)
                [~,Best_ind] = min(obj.lower_CVs);
                value = obj(Best_ind);
            elseif length(obj(index(1)).lower_obj) > 1
                Best = NDSort(obj(index).lower_objs,1);
                Best_ind = index(Best==1);
                value = obj(Best_ind);
            else
                [~,Best] = min(obj(index).lower_objs);
                Best_ind = index(Best);
                value = obj(Best_ind);
            end
        end
        
        function obj = adds(obj,AddProper,type)
        %adds - Get the matrix of additional properties of the population.
        %
        %   A = obj.adds(AddProper) returns the matrix of additional
        %   properties of the population obj. If any individual in obj does
        %   not contain an additional property, assign it a default value
        %   specified in AddProper.
            if nargin>1
                switch type
                    case 'upper'
                        for i = 1 : length(obj)
                            obj(i).upper_add = AddProper(i,:);
                        end
                    case 'lower'
                        for i = 1 : length(obj)
                            obj(i).lower_add = AddProper(i,:);
                        end
                end
            end
        end
        
        function value = upper_adds(obj,dim)
            if ~isempty(obj)
                value = cat(1,obj.upper_add);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
        
        function value = lower_adds(obj,dim)
            if ~isempty(obj)
                value = cat(1,obj.lower_add);
                if nargin > 1
                    if max(dim) <= size(value,2)
                        value = value(:,dim);
                    else
                        error('dim <= %d',size(value,2));
                    end
                end
            else
                value = [];
            end
        end
    end
end