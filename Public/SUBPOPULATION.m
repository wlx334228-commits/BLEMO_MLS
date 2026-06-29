classdef SUBPOPULATION< handle
    properties(SetAccess = private)
        upper_dec; 
        Individuals;
        Best_Individuals;
        FrontNo;
        CrowdDis;
    end
    
    methods
        function obj = SUBPOPULATION(upper_decs,Individuals,Best_Individuals,FrontNo,CrowdDis)
            if nargin >0
                L = size(upper_decs,1);
                if L == 1
                    obj.upper_dec = upper_decs;
                    obj.Individuals = Individuals;
                    obj.Best_Individuals = Best_Individuals;
                    obj.FrontNo = FrontNo;
                    obj.CrowdDis = CrowdDis;
                else
                    obj(1,L) = SUBPOPULATION;
                    for i=1:L
                        obj(i).upper_dec = upper_decs(i,:);
                        obj(i).Individuals = Individuals{i};
                        obj(i).Best_Individuals = Best_Individuals{i};
                        obj.FrontNo = FrontNo{i};
                        obj.CrowdDis = CrowdDis{i};
                    end
                end
            end
        end
        
        function value = upper_decs(obj)
            if ~isempty(obj)
                value = cat(1,obj.upper_dec);
            else
                value = [];
            end
        end
        
        function [value,index] = Population(obj,type)
            value = [];
            if ~isempty(obj)
                if nargin > 1
                    switch type 
                        case 'cell'
                            for i = 1:length(obj)
                                value = cat(1,value,{obj(i).Individuals});
                            end
                            index = [];
                    end
                else
                    N = length(obj);
                    index = [];
                    value = [];
                    for i=1:N
                        value_i = obj(i).Individuals;
                        value = cat(2,value,value_i);
                        index = cat(2,index,repmat(i,1,length(value_i)));
                    end
                end
            end
        end
        
        function [value,index] = Best_Population(obj,type)
            value = [];
            if ~isempty(obj)
                if nargin > 1
                    switch type 
                        case 'cell'
                            for i = 1:length(obj)
                                value = cat(1,value,{obj(i).Best_Individuals});
                            end
                            index = [];
                    end
                else
                    N = length(obj);
                    index = [];
                    value = [];
                    for i=1:N
                        value_i = obj(i).Best_Individuals;
                        value = cat(2,value,value_i);
                        index = cat(2,index,repmat(i,1,length(value_i)));
                    end
                end
            end
        end
        
        function [best,index] = Bilevel_Best(obj)
            if ~isempty(obj)
                N = length(obj);
                index = [];
                Pop = [];
                for i=1:N
                    Pop_i = obj(i).Best_Individuals;
                    Pop = cat(2,Pop,Pop_i);
                    index = cat(2,index,repmat(i,1,length(Pop_i)));
                end
                [best,inds] = Pop.upper_best;
                index = index(inds);
                for i=1:length(best)
                    best(i).label = true;
                end
            else
                best = [];
                index = [];
            end
        end
    end
end