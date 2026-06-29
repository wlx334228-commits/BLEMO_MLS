classdef DATA < handle
    properties(SetAccess = private)
        upper_dec; 
        Archive;
        Best_Archive;
    end
    
    methods
        function obj = DATA(upper_decs,Archive,Best_Archive)
            if nargin >0
                obj(1,length(Archive)) = DATA;
                for i=1:length(Archive)
                    obj(i).upper_dec = upper_decs(i,:);
                    obj(i).Archive = Archive{i};
                    obj(i).Best_Archive = Best_Archive{i};
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
        
        function value = Archives(obj,type)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            value = [];
            if ~isempty(obj)
%                 for i = 1:length(obj)
%                     value = cat(1,value,{obj(i).Archive});
%                 end
                if nargin > 1
                    switch type 
                        case 'cell'
                            for i = 1:length(obj)
                                value = cat(1,value,{obj(i).Archive});
                            end
                    end
                else
                    value = cat(2,obj.Archive);
                end
            end
        end
        
        function [value,index] = Best_Archives(obj,type)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            value = [];
            if ~isempty(obj)
%                 for i = 1:length(obj)
%                     value = cat(1,value,{obj(i).Best_Archive});
%                 end
                if nargin > 1
                    switch type 
                        case 'cell'
                            for i = 1:length(obj)
                                value = cat(1,value,{obj(i).Archive});
                            end
                            index = [];
                    end
                else
                    N = length(obj);
                    index = [];
                    value = [];
                    for i=1:N
                        value_i = obj(i).Best_Archive;
                        value = cat(2,value,value_i);
                        index = cat(2,index,repmat(i,1,length(value_i)));
                    end
                end
            end
        end
        
        function [best,index] = bilevel_best(obj)
            if ~isempty(obj)
                N = length(obj);
                index = [];
                Pop = [];
                for i=1:N
                    Pop_i = obj(i).Best_Archive;
                    Pop = cat(2,Pop,Pop_i);
                    index = cat(2,index,repmat(i,1,length(Pop_i)));
                end
                [best,inds] = Pop.bilevel_best;
                index = index(inds);
                for i=1:length(best)
                    best(i).label = true;
                end
%                 best = DATA(value,obj(index).Archives,obj(index).Best_Archives);
            else
                best = [];
                index = [];
            end
        end
    end
end