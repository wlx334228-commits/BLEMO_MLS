classdef GLOBAL < handle
    properties
        N;                  % Population size
    end
    properties(SetAccess = private)
        problem;            % Problem function
        algorithm;
        maxFEs;
        upper_FEs = 0;                 % Number of evaluated individuals
        lower_FEs = 0; 
        PF = [];
        gen=0;                            % Current generation
        run        = 1;                 % Run number
        runtime    = 0;                 % Runtime
        save       = 0;             	% Number of saved populations
        result     ={};
        CalMetric = {'runtime','upper_FEs','lower_FEs'}
        Fig = struct('ShowFig',false,'Figtype',1,'currentFig',[]);
        parameter  = struct();      	% Parameters of functions specified by users
        outputFcn  = @GLOBAL.Output;  	% Function invoked after each generation
        isparallel = false;
    end
    
    properties(SetAccess = ?PROBLEM)
        M;                              % Number of objectives
        D;                              % Number of decision variables
        upper_domain;                         % Upper domain of each decision variable
        lower_domain;                         % Lower domain of each decision variable
        Parameter;
    end
    
    %%
    methods
        function obj = GLOBAL(varargin)
            obj.GetObj(obj);
            % Initialize the parameters which can be specified by users
            propertyStr = {'algorithm','problem','CalMetric','N','M','D','maxgen','maxFEs','run','save','Fig','isparallel'};
            if nargin > 0
                IsString = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
                [~,Loc]  = ismember(varargin(IsString),cellfun(@(S)['-',S],propertyStr,'UniformOutput',false));
                for i = Loc
                    if i==1
                        obj.(propertyStr{i}) = str2func(varargin{IsString(Loc==i)+1}{1});
                    elseif i==2
                        obj.(propertyStr{i}) = str2func(varargin{IsString(Loc==i)+1}{1});
                        obj.D = varargin{IsString(Loc==i)+1}{2};
                        obj.N = varargin{IsString(Loc==i)+1}{3};
                    elseif i==3
                        obj.(propertyStr{i}) = cat(2,obj.(propertyStr{i}),varargin{IsString(Loc==i)+1});
                    else
                        obj.(propertyStr{i}) = varargin{IsString(Loc==i)+1};
                    end
                end
            end
            % Instantiate a problem object
            obj.problem = obj.problem();
            obj.PF = obj.problem.PF();
            % Add the folders of the algorithm and problem to the top of
            % the search path
            addpath(fileparts(which(class(obj.problem))));
            addpath(fileparts(which(func2str(obj.algorithm))));
        end
        
        function Start(obj)
            folder = fullfile('Data',func2str(obj.algorithm));
            if obj.isparallel && exist(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),obj.run)),'file')
                return;
            end
            if obj.upper_FEs <= 0
                if obj.Fig.ShowFig
                    close('all');
                    obj.Fig.currentFig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_R%d',func2str(obj.algorithm),class(obj.problem),obj.M(1),obj.M(2),obj.D(1),obj.D(2),obj.run),'NumberTitle','off','Position',[0,0,900,450],'Visible','off');
                    movegui(obj.Fig.currentFig,'north');
                end
                
                try
                    tic;
                    warning off;
                    obj.algorithm(obj);
                catch err
                    if strcmp(err.identifier,'GLOBAL:Termination')||obj.isparallel
                        return;
                    else
                        rethrow(err);
                    end
                end
            end
        end
        
        function notermination = NotTermination(obj,Population,notermination)
            
            % Accumulate the runtime
            obj.runtime = obj.runtime + toc;
            
            % Save the last population
            if obj.save>0
                if mod(obj.gen-1,obj.save) == 0
                    index = size(obj.result,1)+1;
                else
                    index = max(1,size(obj.result,1));
                end
            else
                index = 1;
            end
            
            obj.result(index,:).upper_FEs = obj.upper_FEs;
            obj.result(index,:).lower_FEs = obj.lower_FEs;
            obj.result(index,:).Population = Population;
            
            % Invoke obj.outputFcn
            drawnow();
            
            % Detect whether the total number of evaluations has exceeded
            if nargin < 3
                if obj.upper_FEs + obj.lower_FEs >= obj.TotalMaxFEs()
                    notermination = false;
                else
                    notermination = true;
                end
            end
            
            if obj.upper_FEs > 0
                obj.outputFcn(obj,notermination);
            end
            
            obj.gen = obj.gen +1;
            
            assert(notermination,'GLOBAL:Termination','Algorithm has terminated');
            
            tic;
        end

        function value = TotalMaxFEs(obj)
            value = obj.maxFEs;
            if isempty(value)
                value = inf;
            elseif numel(value) > 1
                value = sum(value(:));
            end
        end
        
       %% 初始化
        function Population = Initialization(obj,parameter,N)
            if nargin < 2
                N = obj.N(1);
                domain = [obj.upper_domain,obj.lower_domain];
                PopDec = unifrnd(repmat(domain(1,:),N,1),repmat(domain(2,:),N,1));
                Upper_decs = obj.problem.Decs(PopDec(:,1:size(obj.upper_domain,2)),'upper');
                Lower_decs = obj.problem.Decs(PopDec(:,size(obj.upper_domain,2)+1:end),'lower');
                Population = INDIVIDUAL({Upper_decs,Lower_decs});
            else
                if isa(parameter,'char')
                    switch parameter
                        case 'upper'
                            if nargin ==2
                                N = obj.N(1);
                            end
                            domain = obj.upper_domain;
                            PopDec = unifrnd(repmat(domain(1,:),N,1),repmat(domain(2,:),N,1));
                            PopDec = obj.problem.Decs(PopDec,parameter);
                            Population = INDIVIDUAL(PopDec,'upper');
                        case 'lower'
                            if nargin ==2
                                N = obj.N(2);
                            end
                            domain = obj.lower_domain;
                            PopDec = unifrnd(repmat(domain(1,:),N,1),repmat(domain(2,:),N,1));
                            PopDec = obj.problem.Decs(PopDec,parameter);
                            Population = INDIVIDUAL(PopDec,'lower');
                        case 'bilevel'
                            if nargin ==2
                                N = obj.N(2);
                            end
                            domain = [obj.upper_domain,obj.lower_domain];
                            Dec = unifrnd(repmat(domain(1,:),N,1),repmat(domain(2,:),N,1));
                            upper_Dec = obj.problem.Decs(Dec(:,1:obj.D(1)),'upper');
                            lower_Dec = obj.problem.Decs(Dec(:,1+obj.D(1):end),'lower');
                            Population = INDIVIDUAL({upper_Dec,lower_Dec});
                    end
                elseif isa(parameter,'double')
                    if size(parameter,2)==obj.D(1)
                        upper_dec = parameter;
                        upper_dec = obj.problem.Decs(upper_dec,'upper');
                        if nargin ==2
                            N = obj.N(2);
                        end
                        domain = obj.lower_domain;
                        PopDec = unifrnd(repmat(domain(1,:),N,1),repmat(domain(2,:),N,1));
                        PopDec = obj.problem.Decs(PopDec,'lower');
                        Population = INDIVIDUAL({upper_dec,PopDec});
                    elseif sum(size(parameter))==2
                        N = parameter;
                        domain = [obj.upper_domain,obj.lower_domain];
                        PopDec = unifrnd(repmat(domain(1,:),N,1),repmat(domain(2,:),N,1));
                        Upper_decs = obj.problem.Decs(PopDec(:,1:size(obj.upper_domain,2)),'upper');
                        Lower_decs = obj.problem.Decs(PopDec(:,size(obj.upper_domain,2)+1:end),'lower');
                        Population = INDIVIDUAL({Upper_decs,Lower_decs});
                    end
                end
            end
        end
        %% 获取给定上层变量对应的下层PF
        function lower_PF = get_lower_PF(obj,upper_dec)
            lower_PF = obj.problem.lower_PF(upper_dec);
        end
        
        function SetFEs(obj,FEs)
            obj.upper_FEs = FEs(1);
            obj.lower_FEs = FEs(2);
        end
        
        function value = FEs(obj)
            value = [obj.upper_FEs,obj.lower_FEs];
        end
        
       %% 评估目标函数
        function Population = Evaluate(obj,Population,type,target)
            if ~isempty(Population)
                if nargin <= 3
                    target = {'obj' 'con'};
                    if nargin <= 2
                        type = {'upper' 'lower'};
                    else
                        type = cellstr(type);
                    end
                else
                    target = cellstr(target);
                    type = cellstr(type);
                end
                
                for k = 1:length(type)
                    switch type{k}
                        case 'upper'
                            for l=1:length(target)
                                switch target{l}
                                    case 'obj'
                                        index = find(~Population.upper_CalObjs);
                                        if ~isempty(index)
                                            try 
                                                [objs,adds] = obj.problem.CalObj(Population(index),'upper');
                                                Population(index) = Population(index).update_value(adds,'add');
                                            catch
                                                objs = obj.problem.CalObj(Population(index),'upper');
                                            end
                                            obj.upper_FEs = obj.upper_FEs + length(index);
                                            Population(index) = Population(index).update_value(objs,'upper','obj');
                                        end
                                    case 'con'
                                        index = find(~Population.upper_CalCons);
                                        if ~isempty(index)
                                            cons = obj.problem.CalCon(Population(index),'upper');
                                            Population(index) = Population(index).update_value(cons,'upper','con');
                                        end
                                end
                            end
                        case 'lower'
                            for l=1:length(target)
                                switch target{l}
                                    case 'obj'
                                        index = find(~Population.lower_CalObjs);
                                        if ~isempty(index)
                                            try 
                                                [objs,adds] = obj.problem.CalObj(Population(index),'lower');
                                                Population(index) = Population(index).update_value(adds,'add');
                                            catch
                                                objs = obj.problem.CalObj(Population(index),'lower');
                                            end
                                            obj.lower_FEs = obj.lower_FEs + length(index);
                                            Population(index) = Population(index).update_value(objs,'lower','obj');
                                        end
                                    case 'con'
                                        index = find(~Population.lower_CalCons);
                                        if ~isempty(index)
                                            cons = obj.problem.CalCon(Population(index),'lower');
                                            Population(index) = Population(index).update_value(cons,'lower','con');
                                        end
                                end
                            end
                    end
                end
            end
        end
        
        
        %% 获取参数
        function varargout = ParameterSet(obj,varargin)
        %ParameterSet - Obtain the parameter settings from user.
        %
        %   [p1,p2,...] = obj.ParameterSet(v1,v2,...) returns the values of
        %   p1, p2, ..., where v1, v2, ... are their default values. The
        %   values are specified by the user with the following form:
        %   MOEA(...,'-X_parameter',{p1,p2,...},...), where X is the
        %   function name of the caller.
        %
        %   MOEA(...,'-X_parameter',{[],p2,...},...) indicates that p1 is
        %   not specified by the user, and p1 equals to its default value
        %   v1.
        %
        %   Example:
        %       [p1,p2,p3] = obj.ParameterSet(1,2,3)

            CallStack = dbstack();
            caller    = CallStack(2).file;
            caller    = caller(1:end-2);
            varargout = varargin;
            if isfield(obj.parameter,caller)
                specified = cellfun(@(S)~isempty(S),obj.parameter.(caller));
                varargout(specified) = obj.parameter.(caller)(specified);
            end
        end
        
        %% 输出下层优化
        function lower_Output(obj,Pop,pos)
            LPF = [];
            if nargin == 2
                pos = 4;
            end
            if obj.Fig.ShowFig
                try
                    obj.Fig.currentFig = figure(obj.Fig.currentFig);
                catch
                    obj.Fig.currentFig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_R%d',func2str(obj.algorithm),class(obj.problem),obj.M(1),obj.M(2),obj.D(1),obj.D(2),obj.run),'NumberTitle','off','Position',[0,0,900,450]);
                    movegui(obj.Fig.currentFig,'north');
                end
                if obj.Fig.Figtype == 1
                    current_subplot = subplot(2,3,min(pos,3),'Parent',obj.Fig.currentFig,'FontName','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                elseif obj.Fig.Figtype == 2
                    current_subplot = subplot(2,2,min(pos,2),'Parent',obj.Fig.currentFig,'FontName','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                end
                if obj.M(2)>1
                    set(gcf,'CurrentAxes',current_subplot)
                    if ~isempty(LPF)
                        plot(LPF(:,1),LPF(:,2),'Parent',current_subplot,'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                        set(current_subplot,'NextPlot','add')
                    end
                    lower_feasibles = Pop.lower_feasibles;
                    if isempty(Pop(lower_feasibles))
                        set(current_subplot,'NextPlot','replacechildren');
                        plot(Pop(~lower_feasibles).lower_objs(1),Pop(~lower_feasibles).lower_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    else
                        plot(Pop(lower_feasibles).lower_objs(1),Pop(lower_feasibles).lower_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                        set(current_subplot,'NextPlot','add')
                        plot(Pop(~lower_feasibles).lower_objs(1),Pop(~lower_feasibles).lower_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    end
                    xlabel('$f_1$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    ylabel('$f_2$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    title('Lower Objective Space','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    
                else
                    set(gcf,'CurrentAxes',current_subplot)
                    if ~isempty(LPF)
                        plot(1:length(LPF),LPF(:,1),'Parent',current_subplot,'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                        set(current_subplot,'NextPlot','add')
                    end
                    ii = 1:length(Pop);
                    lower_feasibles = Pop.lower_feasibles;
                    if isempty(Pop(lower_feasibles))
                        set(current_subplot,'NextPlot','replacechildren');
                        plot(ii(~lower_feasibles),Pop(~lower_feasibles).lower_objs(1),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    else
                        plot(ii(lower_feasibles),Pop(lower_feasibles).lower_objs(1),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                        set(current_subplot,'NextPlot','add')
                        plot(ii(~lower_feasibles),Pop(~lower_feasibles).lower_objs(1),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    end
                    xlabel('$f_1$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    ylabel('$f_2$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    title('Lower Objective Space','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                end
                if obj.Fig.Figtype == 1
                    current_subplot = subplot(2,3,min(pos+3,6),'Parent',obj.Fig.currentFig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                elseif obj.Fig.Figtype == 2
                    current_subplot = subplot(2,2,min(pos+2,4),'Parent',obj.Fig.currentFig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                end
                set(gcf,'CurrentAxes',current_subplot)
                plot(transpose(Pop.decs));
                xlim([1 size(Pop.decs,2)]);
                xlabel('Dim','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                ylabel('Value','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                title('Decision Space','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
            end
        end
        
        %% 输出上层优化
        function upper_Output(obj,Pop,pos)
            
            if nargin == 2
                pos = 3;
                name = inputname(2);
            elseif nargin == 3
                name = inputname(2);
            end
            
            if isempty(name)
                name = 'Upper Objective Space';
            end
            
            if obj.Fig.ShowFig
                try
                    obj.Fig.currentFig = figure(obj.Fig.currentFig);
                catch
                    obj.Fig.currentFig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_R%d',func2str(obj.algorithm),class(obj.problem),obj.M(1),obj.M(2),obj.D(1),obj.D(2),obj.run),'NumberTitle','off','Position',[0,0,900,450]);
                    movegui(obj.Fig.currentFig,'north');
                end
                
                if obj.Fig.Figtype == 1
                    current_subplot = subplot(2,3,min(pos,3),'Parent',obj.Fig.currentFig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                elseif obj.Fig.Figtype == 2
                    current_subplot = subplot(2,2,min(pos,2),'Parent',obj.Fig.currentFig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                end
                if obj.M(2)>1
                    set(gcf,'CurrentAxes',current_subplot)
                    if ~isempty(obj.PF)
                        plot(obj.PF(:,1),obj.PF(:,2),'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                        set(current_subplot,'NextPlot','add')
                    end
                    
                    statuses = Pop.upper_CalObjs & Pop.lower_CalObjs;
                    Pop_temp = Pop(statuses);
                    feasibles = Pop_temp.upper_feasibles & Pop_temp.lower_feasibles;
                    if any(feasibles)
                        plot(Pop_temp(feasibles).upper_objs(1),Pop_temp(feasibles).upper_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                        set(current_subplot,'NextPlot','add');
                    end
                    
                    infeasibles = (Pop.upper_CalObjs & ~Pop.upper_feasibles)|(Pop.lower_CalObjs & ~Pop.lower_feasibles);
                    if any(infeasibles)
                        plot(Pop(infeasibles).upper_objs(1),Pop(infeasibles).upper_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                        set(current_subplot,'NextPlot','add');
                    end
                    check_infeasibles = (Pop.upper_CalObjs & Pop.upper_feasibles) & ~Pop.lower_CalObjs;
                    plot(Pop(check_infeasibles).upper_objs(1),Pop(check_infeasibles).upper_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','g');
                    
                    xlabel('$F_1$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    ylabel('$F_2$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    title(name,'Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                else
                    set(gcf,'CurrentAxes',current_subplot)
                    if ~isempty(obj.PF)
                        plot(1:length(Pop),repmat(obj.PF,1,length(Pop)),'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                        set(current_subplot,'NextPlot','add')
                    end
                    ii = 1:length(Pop);
                    statuses = Pop.upper_CalObjs & Pop.lower_CalObjs;
                    Pop_temp = Pop(statuses);
                    ii_statuses = ii(statuses);
                    feasibles = Pop_temp.upper_feasibles & Pop_temp.lower_feasibles;
                    plot(ii_statuses(feasibles),Pop_temp(feasibles).upper_objs(1),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                    
                    infeasibles = (Pop.upper_CalObjs & ~Pop.upper_feasibles)|(Pop.lower_CalObjs & ~Pop.lower_feasibles);
                    plot(ii(infeasibles),Pop(infeasibles).upper_objs(1),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    
                    check_infeasibles = (Pop.upper_CalObjs & Pop.upper_feasibles) & ~Pop.lower_CalObjs;
                    plot(ii(check_infeasibles),Pop(check_infeasibles).upper_objs(1),'MarkerSize',6,'Marker','o','LineStyle','none','color','g');
                    
                    xlabel('$F_1$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    ylabel('$F_2$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                    title(name,'Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                end
                
                if obj.Fig.Figtype == 1
                    current_subplot = subplot(2,3,min(pos,3)+3,'Parent',obj.Fig.currentFig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                elseif obj.Fig.Figtype == 2
                    current_subplot = subplot(2,2,min(pos,3)+2,'Parent',obj.Fig.currentFig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                end
                set(gcf,'CurrentAxes',current_subplot)
                
                plot(transpose(Pop.decs),'Parent',current_subplot);
                xlim([1 size(Pop.decs,2)]);
                xlabel('Dim','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                ylabel('Value','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                title('Decision Space','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
            end
        end 
        
    end
    
    methods(Access = private, Static)
        %% Display or save the result after the algorithm is terminated
        function Output(obj,notermination)
            if ~notermination
                clc;
            end
            if obj.gen > 0
                fprintf('R%d-G%d: %s on %s, %d+%d objectives %d+%d variables, FEs=[%d,%d], %.2fs passed...\n',...
                    obj.run,obj.gen,func2str(obj.algorithm),class(obj.problem),obj.M(1),obj.M(2),obj.D(1),obj.D(2),obj.upper_FEs,obj.lower_FEs,obj.runtime);
            end
            
            if ~isempty(obj.result(end).Population)
                obj.upper_Output(obj.result(end).Population,1);
            end
            
            if ~notermination
                if obj.Fig.ShowFig
                    folder = fullfile('Figures',func2str(obj.algorithm));
                    [~,~]  = mkdir(folder);
                    saveas(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.fig',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),obj.run)));
                    %                     print(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_R%d',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),obj.run)),'-depsc');
                    %                     saveas(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.svg',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),obj.run)));
                    close('all');
                else
                    GLOBAL.SaveFig(obj);
                end
            end
            
            if ~notermination
                folder = fullfile('Data',func2str(obj.algorithm));
                [~,~]  = mkdir(folder);
                Data = obj;
                Result.runtime = Data.runtime;
                Result.upper_FEs = Data.upper_FEs;
                Result.lower_FEs = Data.lower_FEs;
                for i=1:length(Data.result)
                    Population_i = Data.result(i).Population;
                    if ~isempty(Population_i)
                        Feasible     = Population_i.bilevel_feasibles;
                        Population_i = Population_i(Feasible);
                    end
                    if ~isempty(Population_i) && ~isempty(Data.PF) && size(Data.PF,2)>1
                        Population_i(any(Population_i.upper_objs>repmat(max(Data.PF,[],1),length(Population_i),1),2))=[];
                        if ~isempty(Population_i)
                            [~,ID] = min(pdist2(Population_i.upper_objs,Data.PF),[],2);
                            Index = all(Population_i.upper_objs<Data.PF(ID,:),2);
                            Population_i = Population_i(~Index);
                        end
                    end
                    
                    if length(obj.CalMetric) > 3
                        for j=4:length(obj.CalMetric)
                            metricName = obj.CalMetric{j};
                            if strcmp(metricName,'LGD')
                                try
                                    Feasible     = Data.result(i).Population.lower_feasibles;
                                    Data.result(i).(metricName) = GLOBAL.Metric(str2func(metricName),Data.result(i).Population(Feasible).upper_best,Data);
                                catch
                                    Data.result(i).(metricName) = NaN;
                                end
                            elseif strcmp(metricName,'HV')
                                try
                                    Data.result(i).(metricName) = GLOBAL.Metric(str2func(metricName),Data.result(i).Population,Data);
                                catch
                                    Data.result(i).(metricName) = NaN;
                                end
                            else
                                try
                                    Data.result(i).(metricName) = GLOBAL.Metric(str2func(metricName),Population_i.upper_best,Data);
                                catch
                                    Data.result(i).(metricName) = NaN;
                                end
                            end
                            if i==length(Data.result)
                                Result.(metricName) = Data.result(i).(metricName);
                            end
                        end
                    end
                end
                
                save(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),obj.run)),'Data','Result');
                
            end
            
        end
        
        function value = Metric(metric,Population,Global)
            % Calculate the metric value of the population
            try
                value = metric(Population,Global);
            catch
                try 
                    value = metric(Population.upper_obj,Global.PF);
                catch 
                    value = NaN;
                end
            end
        end
        
    end
    
    methods(Static)
        %% Get the current GLOBAL object
        function obj = GetObj(obj)
        %GetObj - Get the current GLOBAL object.
        %
        %   Global = GLOBAL.GetObj() returns the current GLOBAL object.
        %
        %   Example:
        %       Global = GLOBAL.GetObj()
        
            persistent Global;
            if nargin > 0
                Global = obj;
            else
                obj = Global;
            end
        end
        
        function extract_data(Problems,Algorithms,maxFEs,Runs)
             for j=1:size(Problems,1)
                    Problems_name = Problems{j,1};
                    ProblemsD = Problems{j,2};
                    for i=1:length(Algorithms)
                        folder = fullfile('Data',Algorithms{i});
                        [~,~]  = mkdir(folder);
                        for k=1:Runs
                            try
                                load(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',Algorithms{i},Problems_name,ProblemsD(1),ProblemsD(2),k)));
                                if ~strcmp(Algorithms{i},'BLEADPL')
                                    inds1 = [Data.result.upper_FEs]>= maxFEs{j}(1);
                                    inds1(sum(~inds1)+1) = false;
                                    Data.result(inds1)=[];
                                else
                                    inds2 = [Data.result.lower_FEs]>= maxFEs{j}(2);
                                    inds2(sum(~inds2)+1) = false;
                                    Data.result(inds2)=[];
                                end
                                Result.upper_FEs = Data.result(end).upper_FEs;
                                Result.lower_FEs = Data.result(end).lower_FEs;
                                for p=4:length(Data.CalMetric)
                                    metricName = Data.CalMetric{p};
                                    Result.(metricName) = Data.result(end).(metricName);
                                end
                                save(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',func2str(Data.algorithm),class(Data.problem),Data.D(1),Data.D(2),Data.run)),'Data','Result');
                            catch
                                continue
                            end
                        end
                    end
             end
        end
        
        function Data_dealing(Problems,Algorithms,Run,Metrics,type,RecalM,ReCPF,ReCMetrics)
            metric_names = {'runtime','upper_FEs','lower_FEs','FEs'};
            if nargin > 3
                metric_names = cat(2,metric_names,Metrics);
            end
            
            if nargin < 6
                RecalM = false;
                ReCPF = false;
            end
            
            if RecalM
                if nargin < 7
                    ReCPF = false;
                end
                if nargin < 8
                    ReCMetrics = Metrics;
                end
            end
            
            if nargin<5
                type = 'std';
            end
            
            if ReCPF
                PF_new = {};
                for j=1:size(Problems,1)
                    Problems_name = Problems{j,1};
                    ProblemsD = Problems{j,2};
                    PF_j = [];
                    for i=1:length(Algorithms)
                        folder = fullfile('Data',Algorithms{i});
                        [~,~]  = mkdir(folder);
                        for k=1:Run
                            try
                                load(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',Algorithms{i},Problems_name,ProblemsD(1),ProblemsD(2),k)));
                                PF_j = cat(2,PF_j,Data.result(end).Population);
                            catch
                                continue
                            end
                        end
                    end
                    PF_j = [Data.problem.PF();PF_j.bilevel_best.upper_objs];
                    PF_new = cat(2,PF_new,{PF_j});
                end
            end
            
            for i=1:length(Algorithms)
                
                Title = cell(2,2+3*length(metric_names));
                Title(1,4:3:3*length(metric_names)+1) = metric_names;
                Title(2,:) = [{'Algorithms','Problems'},repmat({'Best result','Median result','Worst result'},1,length(metric_names))];
                
                folder = fullfile('Data',Algorithms{i});
                [~,~]  = mkdir(folder);
                time = clock;
                File_Name = fullfile(folder,sprintf('Result_of_%s_%s.xlsx',Algorithms{i},datestr(time,'mm_dd_HH_MM')));
                
                writecell([{'Algorithms','Problems'},metric_names],File_Name,'Sheet','Mean_Result','WriteMode','overwritesheet');
                writecell(Title,File_Name,'Sheet','Best_Worst_Median_Result','WriteMode','overwritesheet');
                
                write_result.mean_result ={};
                write_result.Best_Worst_Median_Result = {};
                
                for j=1:size(Problems,1)
                    
                    Problems_name = Problems{j,1};
                    ProblemsD = Problems{j,2};
                    
                    for k=1:length(metric_names)
                        result.(metric_names{k}) = [];
                    end
                    
                    Global = [];
%                     Population = [];
                    
                    for k=1:Run
                        try
                            load(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',Algorithms{i},Problems_name,ProblemsD(1),ProblemsD(2),k)));
                            
                            %                         if size(Data.PF,1)<=1 && Data.M(1)>1
                            %                             RecalM = true;
                            %                             Population = cat(2,Population,Data.result(end).Population);
                            %                         end
                            
                            Result.FEs = Result.upper_FEs+Result.lower_FEs;
                            
                            if ~RecalM
                                for l=1:length(metric_names)
                                    if ~isempty(Result.(metric_names{l}))
                                        result.(metric_names{l})=cat(2,result.(metric_names{l}),Result.(metric_names{l}));
                                    else
                                        result.(metric_names{l})=cat(2,result.(metric_names{l}),NaN);
                                    end
                                end
                            end
                            Global = cat(2,Global,Data);
                        catch
                            continue
                        end
                    end
                    
                    if RecalM

%                         if ~isempty(Population)
%                             PF = [Data.problem.PF();Population.upper_best.upper_objs];
%                         else
%                             try 
%                                 PF = Data.problem.PF;
%                             catch
%                                 [~,~,PF,~]=getOptimalSolutionSMD(Data.D(1),Data.D(2),func2str(Data.problem));
%                             end
%                         end

                        for k=1:length(Global)
                            if ReCPF
                                Global(k).PF = PF_new{j};
                            else
                                Global(k).PF = Data.problem.PF();
                            end
                            
                            [Global(k),Result] = GLOBAL.CalMetrics(Global(k),Metrics,ReCMetrics);
                            Result.FEs = Result.upper_FEs+Result.lower_FEs;

                            for l=1:length(metric_names)
                                if ~isempty(Result.(metric_names{l}))
                                    result.(metric_names{l})=cat(2,result.(metric_names{l}),Result.(metric_names{l}));
                                else
                                    result.(metric_names{l})=cat(2,result.(metric_names{l}),NaN);
                                end
                            end
                        end
                    end
                    
                    if length(Algorithms)>1 && i==1
                        ranksum_result.(Problems{j,1}) = result;
                    end
                    
%                     GLOBAL.SaveFig(Global,'all');
                    
%                     M = transpose(cell2mat(struct2cell(result)));
                    
                    Best_Worst_Median_Result_j = {Algorithms{i},Problems_name};
                    
                    for k=1:length(metric_names)
                        result_k = result.(metric_names{k});
                        ID = find(~isnan(result_k));
                        if ~isempty(ID)
                            [~,rank]= sort(result_k(ID),'ascend');
                            if strcmp(metric_names{k},'UHV') || strcmp(metric_names{k},'LHV')
                                Best_ind = ID(rank(end));
                                Worst_ind = ID(rank(1));
                            else
                                Best_ind = ID(rank(1));
                                Worst_ind = ID(rank(end));
                            end
                            Median_ind = ID(rank(ceil(length(ID)/2)));
                        else
                            IDrand = randperm(length(result_k),3);
                            Best_ind = IDrand(1);
                            Worst_ind = IDrand(2);
                            Median_ind = IDrand(3);
                        end
                        
                        Best_Worst_Median_Result_j = cat(2,Best_Worst_Median_Result_j,{result_k(Best_ind),result_k(Median_ind),result_k(Worst_ind)});

                        index = cellfun(@(x) strcmp(metric_names{k},x),Metrics);
                        if any(index)
                            GLOBAL.SaveFig(Global(Best_ind),sprintf('Best_%s',metric_names{k}));
                            GLOBAL.SaveFig(Global(Worst_ind),sprintf('Worst_%s',metric_names{k}));
                            GLOBAL.SaveFig(Global(Median_ind),sprintf('Median_%s',metric_names{k}));
%                             
%                             GLOBAL.SaveFig(Global(Best_ind),sprintf('Best_%s',metric_names{k}),[]);
%                             GLOBAL.SaveFig(Global(Worst_ind),sprintf('Worst_%s',metric_names{k}),[]);
%                             GLOBAL.SaveFig(Global(Median_ind),sprintf('Median_%s',metric_names{k}),[]);
                        end
                        
                        switch type
                            case 'std'
                                Vars = std(result_k(ID));
                            case 'iqr'
                                Vars = iqr(result_k(ID));
                        end
                        
                        if ~isempty(ID)
                            
                            if length(Algorithms)>1 && i>1
                                ranksum_result_j = ranksum_result.(Problems{j,1});
                                ranksum_result_jk = ranksum_result_j.(metric_names{k});
                            end
                            
                            if length(Algorithms)>1 && i>1 && k>4
                                if any(~isnan(ranksum_result_jk))
                                    [~,p] = ranksum(result_k(ID),ranksum_result_jk);
                                else
                                    p = false;
                                end
                                if p
                                    if strcmp(metric_names{k},'UHV') || strcmp(metric_names{k},'LHV')
                                        if mean(ranksum_result_jk(~isnan(ranksum_result_jk))) < mean(result_k(ID))
                                            mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)%s',mean(result_k(ID)),Vars,'+');
                                            %                                             mean_result.(metric_names{k}) =  sprintf('%0.3d%s',mean(result_k(ID)),'+');
                                        else
                                            mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)%s',mean(result_k(ID)),Vars,'-');
                                            %                                             mean_result.(metric_names{k}) =  sprintf('%0.3d%s',mean(result_k(ID)),'-');
                                        end
                                    else
                                        if mean(ranksum_result_jk(~isnan(ranksum_result_jk))) > mean(result_k(ID))
                                            mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)%s',mean(result_k(ID)),Vars,'+');
                                            %                                             mean_result.(metric_names{k}) =  sprintf('%0.3d%s',mean(result_k(ID)),'+');
                                        else
                                            mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)%s',mean(result_k(ID)),Vars,'-');
                                            %                                             mean_result.(metric_names{k}) =  sprintf('%0.3d%s',mean(result_k(ID)),'-');
                                        end
                                    end
                                else
                                    if all(isnan(ranksum_result_jk))
                                        mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)%s',mean(result_k(ID)),Vars,'-');
                                    else
                                        mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)%s',mean(result_k(ID)),Vars,'=');
                                    end
                                    %                                     mean_result.(metric_names{k}) =  sprintf('%0.3d%s',mean(result_k(ID)),'=');
                                end
                            else
                                if i>1 && k>1
%                                     mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)[\\textbf{%0.2f}]',mean(result_k(ID)),Vars,mean(result_k(ID))/mean(ranksum_result_jk(~isnan(ranksum_result_jk))));
                                    mean_result.(metric_names{k}) =  sprintf('%0.2e[\\textbf{%0.2f}]',mean(result_k(ID)),mean(result_k(ID))/mean(ranksum_result_jk(~isnan(ranksum_result_jk))));
                                elseif i==1 && k>4
                                    mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)',mean(result_k(ID)),Vars);
                                else
                                    mean_result.(metric_names{k}) = mean(result_k(ID));
                                end
                            end
                        else
                            if length(Algorithms)>1 && i>1 && k>4
                                mean_result.(metric_names{k}) = 'NaN(NaN)-';
%                                 mean_result.(metric_names{k}) = 'NaN+';
                            else
                                if k>4
                                    mean_result.(metric_names{k}) = 'NaN(NaN)';
                                    mean_result.(metric_names{k}) = 'NaN';
                                else
                                    mean_result.(metric_names{k}) = 'NaN';
                                end
                            end
                        end
                        
                    end
                    
                    write_result.Best_Worst_Median_Result = cat(1,write_result.Best_Worst_Median_Result,Best_Worst_Median_Result_j);
                    write_result.mean_result = cat(1,write_result.mean_result,cat(2,{Algorithms{i},Problems_name},transpose(struct2cell(mean_result))));
                    
                end
                
                fields = fieldnames(write_result);
                for ii=1:length(fields)
                    writecell(write_result.(fields{ii}),File_Name,'Sheet',fields{ii},'WriteMode','append');
                end
            end
        end
        
        function PlotPFs(Problems,Algorithms,Run,metric_names,legends)
            
            LGDs = cat(2,{'BPF'},legends(2:end),legends(1));
            Algorithms = cat(1,Algorithms(2:end),Algorithms(1));
            
            for j=1:size(Problems,1)
                
                Problems_name = Problems{j,1};
                ProblemsD = Problems{j,2};
                
                Plots_j = [];
                
                close('all');
                movegui(figure('Name',sprintf('The Plot on %s',Problems_name),'Position',[463,244,657,571],'Visible','on'),'north');
                marks ={'*','p','h','s','d','o'};
                colors = {[0.2, 0.8, 0.8], [0, 0.45, 0.74],[0.93, 0.69, 0.13],[0.49, 0.18, 0.56],[0.47, 0.67, 0.19],'r','#4DBEEE',...
                    [0.728683436587994,0.887284680027887,0.0558483212688232],[0.364401368261552,0.171481076831545,0.795361567918927],[0.613545882666124,0.988061285755543,0.219900778839496],[0.241285612653001,0.622924050374514,0.522927762157550]};
                
                for i=1:length(Algorithms)
                    
                    folder = fullfile('Data',Algorithms{i});
                    [~,~]  = mkdir(folder);
                    
                    for k=1:length(metric_names)
                        result.(metric_names{k}) = [];
                    end
                    
                    Global = [];
%                     Population = [];
                    
                    for k=1:Run
                        try
                            load(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',Algorithms{i},Problems_name,ProblemsD(1),ProblemsD(2),k)));
                            for l=1:length(metric_names)
                                if ~isempty(Result.(metric_names{l}))
                                    result.(metric_names{l})=cat(2,result.(metric_names{l}),Result.(metric_names{l}));
                                else
                                    result.(metric_names{l})=cat(2,result.(metric_names{l}),NaN);
                                end
                            end
                            Global = cat(2,Global,Data);
                        catch
                            continue
                        end
                    end
                    
                    for k=1:length(metric_names)
                        result_k = result.(metric_names{k});
                        ID = find(~isnan(result_k));
                        if ~isempty(ID)
                            [~,rank]= sort(result_k(ID),'ascend');
                            if strcmp(metric_names{k},'UHV') || strcmp(metric_names{k},'LHV')
                                Best_ind = ID(rank(end));
                                Worst_ind = ID(rank(1));
                            else
                                Best_ind = ID(rank(1));
                                Worst_ind = ID(rank(end));
                            end
                            Median_ind = ID(rank(ceil(length(ID)/2)));
                        else
                            IDrand = randperm(length(result_k),3);
                            Best_ind = IDrand(1);
                            Worst_ind = IDrand(2);
                            Median_ind = IDrand(3);
                        end
                        
                        Plots_j = cat(2,Plots_j,Global(Best_ind));
                        
%                         index = cellfun(@(x) strcmp(metric_names{k},x),Metrics);
%                         if any(index)
%                             GLOBAL.SaveFig(Global(Best_ind),sprintf('Best_%s',metric_names{k}));
%                             GLOBAL.SaveFig(Global(Worst_ind),sprintf('Worst_%s',metric_names{k}));
%                             GLOBAL.SaveFig(Global(Median_ind),sprintf('Median_%s',metric_names{k}));
% %                             
% %                             GLOBAL.SaveFig(Global(Best_ind),sprintf('Best_%s',metric_names{k}),[]);
% %                             GLOBAL.SaveFig(Global(Worst_ind),sprintf('Worst_%s',metric_names{k}),[]);
% %                             GLOBAL.SaveFig(Global(Median_ind),sprintf('Median_%s',metric_names{k}),[]);
%                         end 
                    end
                    
                end
                
                LGDs_temp = LGDs;
                k_temp = [];
                
                PF = Plots_j(1).PF;
                plot(PF(:,1),PF(:,2),'MarkerSize',6,'Marker','.','MarkerSize',10,'LineStyle','none','color',[0,0,0]);
                hold on;
                
                xlimits = [min(PF(:,1)),1+max(PF(:,1))];
                ylimits = [min(PF(:,2)),1+max(PF(:,2))];
                
%                 ax1 = gca;
%                 ylimits = get(ax1,'YLim');
%                 set(ax1,'YLim',max([1.3*ylimits;ylimits]),'LineWidth',2,'FontSize',17,'box','on');
                
                for k=1:length(Plots_j)
                    Pop = Plots_j(k).result(end).Population;
                    if ~isempty(Pop)
                        Feasible     = Pop.bilevel_feasibles;
                        Pop = Pop(Feasible);
                    end
                    if ~isempty(Pop) 
                        [~,ID] = min(pdist2(Pop.upper_objs,PF),[],2);
                        Index = all(Pop.upper_objs<PF(ID,:),2);
                        Pop = Pop(~Index);
                    end
                    Pop_Objs = Pop.upper_objs;
                    statuses = Pop.upper_CalObjs & Pop.lower_CalObjs;
                    Pop_Objs_temp = Pop_Objs(statuses,:);
                    feasibles = Pop(statuses).upper_feasibles & Pop(statuses).lower_feasibles;
                    if any(feasibles)
                        plot(Pop_Objs_temp(feasibles,1),Pop_Objs_temp(feasibles,2),'MarkerSize',10,'Marker',marks{k},'LineStyle','none','color',colors{k});
                    else
                        k_temp = cat(2,k_temp,k+1);
                    end
                    
                end
                
                LGDs_temp(k_temp) = [];
                set(gca,'XLim',xlimits,'YLim',ylimits,'LineWidth',2,'FontSize',20,'box','on');
                xlabel('$F_1$','Interpreter','latex','Fontname','Times New Roman','FontSize',20);
                ylabel('$F_2$','Interpreter','latex','Fontname','Times New Roman','FontSize',20);
                title('Upper Objective Space','Interpreter','latex','Fontname','Times New Roman','FontSize',20);
                legend(LGDs_temp{:},'NumColumns',2);
                
                folder = fullfile('Figures','PFPlots');
                [~,~]  = mkdir(folder);
                
                print(gcf,fullfile(folder,Problems_name),'-depsc');
                
                close all;
            end
            
        end
        
        function Convergence_plot(Problems,Algorithms,Metric,Runs,legends)
            close all;
            for j=1:size(Problems,1)
                Problems_name = Problems{j,1};
                ProblemsD = Problems{j,2};
%                 figure('Name',sprintf('The Convergence Plot on %s',Problems_name),'Position',[520,326,880,652],'Visible','on');
                movegui(figure('Name',sprintf('The Convergence Plot on %s',Problems_name),'Position',[463,244,657,571],'Visible','on'),'north');
                marks ={'o','d','s','*','p','h'};
                colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#A2142F','#4DBEEE',...
                    [0.728683436587994,0.887284680027887,0.0558483212688232],[0.364401368261552,0.171481076831545,0.795361567918927],[0.613545882666124,0.988061285755543,0.219900778839496],[0.241285612653001,0.622924050374514,0.522927762157550]};
%                 legends = {'BL-PIEA','SMS-EMOA','H-BLEMO','N-BLEMO','BLEA-DPL'};
%                 legends = {'BL-PIEA','Variant1','Variant2'};
                Picture = [];
                
                for i=1:length(Algorithms)
                    Picture_i = [];
                    for LM = 1:length(Metric)
                        Picture_i.(Metric{LM}) = [];
                        UMetric.(Metric{LM}) = [];
                        FE.(Metric{LM}) = [];
                    end
                    
                    folder = fullfile('Data',Algorithms{i});
                    Global = [];
                    
%                     for k=1:Runs
%                         load(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',Algorithms{i},Problems_name,ProblemsD(1),ProblemsD(2),k)));
%                         Global = cat(2,Global,Data);
%                         for LM = 1:length(Metric)
%                             Metric_k = [Data.result.(Metric{LM})];
%                             FE_k = [Data.result.upper_FEs]+[Data.result.lower_FEs];
%                             try
%                                 UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),Metric_k);
%                                 FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),FE_k);
%                             catch
%                                 Y = size(UMetric.(Metric{LM}),2);
%                                 if length(Metric_k)<Y
%                                     UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),[Metric_k,zeros(1,Y-length(Metric_k))]);
%                                     FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),[FE_k,zeros(1,Y-length(FE_k))]);
%                                 else
%                                     UMetric.(Metric{LM})(:,Y+1:length(Metric_k)) = 0;
%                                     FE.(Metric{LM})(:,Y+1:length(FE_k)) = 0;
%                                     UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),Metric_k);
%                                     FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),FE_k);
%                                 end
%                             end
%                         end
%                     end
                    
                    for k=1:Runs
                        try 
                            load(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',Algorithms{i},Problems_name,ProblemsD(1),ProblemsD(2),k)));
                            Global = cat(2,Global,Data);
                            for LM = 1:length(Metric)
                                Metric_k = [Data.result.(Metric{LM})];
                                FE_k = [Data.result.upper_FEs]+[Data.result.lower_FEs];
                                try
                                    UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),Metric_k);
                                    FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),FE_k);
                                catch
                                    Y = size(UMetric.(Metric{LM}),2);
                                    if length(Metric_k)<Y
                                        UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),[zeros(1,Y-length(Metric_k)),Metric_k]);
                                        FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),[zeros(1,Y-length(FE_k)),FE_k]);
                                    else
                                        UMetric.(Metric{LM}) = cat(2,zeros(size(UMetric.(Metric{LM}),1),length(Metric_k)-Y),UMetric.(Metric{LM}));
                                        FE.(Metric{LM}) = cat(2,zeros(size(FE.(Metric{LM}),1),length(FE_k)-Y),FE.(Metric{LM}));
                                        UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),Metric_k);
                                        FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),FE_k);
                                    end
                                end
                            end
                        catch
                            continue
                        end
                    end
                    
%                     index = all(FE(:,2:end),1);
%                     FEs = [0,mean(FE(:,index),1)];
%                     [~,rank] =sort(FEs,'ascend');
%                     UHV(isnan(UHV))=0;
%                     HVs = [0,mean(UHV(:,index),1)];

                    for LM = 1:length(Metric)
                        
%                         I1 = any(~isnan(UMetric.(Metric{LM})),2);
%                         I2 = any(~isnan(UMetric.(Metric{LM})),1);
                        
                        UMetric_LM = UMetric.(Metric{LM});
                        FEs = FE.(Metric{LM});
                        
                        cFEs = mean(FEs(all(FEs>0,2),:),1);
                        for ic = 1:size(FEs,1)
                            
                            FEs_ic = FEs(ic,FEs(ic,:)>0);
                            UMetric_LM_ic = UMetric_LM(ic,FEs(ic,:)>0);
                            
                            Dis = pdist2(FEs_ic',cFEs');
                            index = zeros(1,length(FEs_ic));
                            index(1)=1;
                            index(end) = length(cFEs);
                            Dis(:,[1 length(cFEs)]) = inf;
                            for ik = 2:length(FEs_ic)-1
                                [~,ind]=min(Dis(ik,:));
                                if ind==1
                                    ind = length(cFEs);
                                end
                                index(ik) = ind;
                                Dis(:,1:ind) = inf;
                            end
                            FEs(ic,:)=0;
                            FEs(ic,index) = FEs_ic;
                            
                            UMetric_LM(ic,:) = NaN;
                            UMetric_LM(ic,index) = UMetric_LM_ic;
                        end
                        
                        FEs = sum(FEs,1)./sum(FEs>0,1);
                        
                        Picture_i.(Metric{LM}).maxFEs = max(FEs);
                        
                        DV = diff(FEs);
                        index = DV<0 & abs(DV)>0.01;
                        while any(index)
                            index = find(index);
                            index = index(index<length(FEs));
                            UMetric_LM(:,index)=[];
                            FEs(index)=[];
                            DV = diff(FEs);
                            index = DV<0 & abs(DV)>0.01;
                        end
                        
                        II =  ~isnan(UMetric_LM);
                        ii = II(:,end)==0;
                        
                        if ~all(ii)
                            UMetric_LM(ii,:)=[];
                            UMetric_LM(isnan(UMetric_LM))=0;
                            UMetric_LM = sum(UMetric_LM,1)./max(sum(II(~ii,:),1),eps);
                        else
                            UMetric_LM = zeros(1,length(FEs));
                        end
                        
%                         II = sum(FEs>0,1)==size(FEs,1);
%                         UMetric_LM = UMetric_LM(:,II);
%                         FEs = mean(FEs(:,II),1);
%                         
%                         II =  ~isnan(UMetric_LM);
%                         I1 = any(II,2)&II(:,end);
%                         I2 = any(II,1);
%                         sumN = sum(II(I1,I2),1);
%                         
%                         if any(I1)
%                             FEs = FEs(:,I2);
%                             UMetric_LM = UMetric_LM(I1,I2);
%                             UMetric_LM(isnan(UMetric_LM))=0;
%                             UMetric_LM = sum(UMetric_LM,1)./(sumN+eps);
%                             
%                             [~,rank] =sort(FEs,'ascend');
%                             FEs = FEs(rank);
%                             UMetric_LM = UMetric_LM(rank);
%                             
%                         else
%                             UMetric_LM = zeros(1,length(FEs));
%                         end
                        
%                         Picture_i.(Metric{LM}).maxFEs = max(FEs);
                        
                        if strcmp(Metric{LM},'UHV')
                            DV = diff(UMetric_LM);
                            index = DV<0 & abs(DV)>0.01;
                            while any(index)
                                index = find(index);
                                index = index(index<length(UMetric_LM));
                                UMetric_LM(index)=[];
                                FEs(index)=[];
                                DV = diff(UMetric_LM);
                                index = DV<0 & abs(DV)>0.01;
                            end
                        else
                            DV = diff(UMetric_LM);
                            index = DV>0 & abs(DV)>0.01;
                            while any(index)
                                index = find(index);
                                UMetric_LM(index)=[];
                                FEs(index)=[];
                                DV = diff(UMetric_LM);
                                index = DV>0 & abs(DV)>0.01;
                            end
                        end
                        
                        Picture_i.(Metric{LM}).FEs = [0,FEs];
                        
%                         if strcmp(Metric{LM},'UHV')
                            Picture_i.(Metric{LM}).value = [0,UMetric_LM];
%                         else
%                             UMetric_LM(UMetric_LM==0) = 1.1;
%                             Picture_i.(Metric{LM}).value = [1.1,UMetric_LM];
%                         end
                    end
                    
                    Picture = cat(1,Picture,Picture_i);
                    
                end
                
                ind = 1;
                
                for LM = 1:length(Metric)
                    Picture_i = [Picture.(Metric{LM})];
                    if strcmp(Metric{LM},'UIGD')
                        max_IGD = max([Picture_i.value]);
                    end
%                     switch Problems_name(1:2)
%                         case 'TP'
%                             maxFEs = median([Picture_i.maxFEs]);
%                         otherwise
%                             maxFEs = median([Picture_i.maxFEs]);
%                     end
                    maxFEs = max([Picture_i.maxFEs]);
                    marks_k = randsample(marks,length(Algorithms));
                    colors_k = randsample(colors,length(Algorithms));
                    
                    for ii=1:length(Picture_i)
                        Picture_ii = Picture_i(ii);
                        if strcmp(Metric{LM},'UIGD')
                            Picture_ii.value(1) = 1.1*max_IGD;
                            Picture_ii.value(Picture_ii.value==0) = max_IGD;
                        end
                        if max(Picture_ii.FEs)<maxFEs
                            tempFEs = max(Picture_ii.FEs):max(Picture_ii.FEs)/length(Picture_ii.FEs):maxFEs;
                            if LM == 1
                                line([Picture_ii.FEs,tempFEs],[Picture_ii.value,repmat(Picture_ii.value(end),1,length(tempFEs))],'Marker',marks_k{ii},'Color',colors{mod(ind,length(Algorithms))+1},'LineWidth',1.25,'MarkerSize',7);
                            else
                                line([Picture_ii.FEs,tempFEs],[Picture_ii.value,repmat(Picture_ii.value(end),1,length(tempFEs))],'Marker',marks_k{ii},'Color',colors{mod(ind,length(Algorithms))+1},'LineStyle','-','LineWidth',1.25,'MarkerSize',7,'Parent',ax2);
                            end
                        else
                            inds = find(Picture_ii.FEs<=maxFEs, 1, 'last' );
                            if inds<length(Picture_ii.FEs)
                                if LM == 1
                                    line([Picture_ii.FEs(1:inds),maxFEs],Picture_ii.value(1:inds+1),'Marker',marks_k{ii},'Color',colors{mod(ind,length(Algorithms))+1},'LineWidth',1.25,'MarkerSize',7);
                                else
                                    line([Picture_ii.FEs(1:inds),maxFEs],Picture_ii.value(1:inds+1),'Marker',marks_k{ii},'Color',colors{mod(ind,length(Algorithms))+1},'LineStyle','-','LineWidth',1.25,'MarkerSize',7,'Parent',ax2);
                                end
                            else
                                if LM == 1
                                    line(Picture_ii.FEs(1:inds),Picture_ii.value(1:inds),'Marker',marks_k{ii},'Color',colors{mod(ind,length(Algorithms))+1},'LineWidth',1.25,'MarkerSize',7);
                                else
                                    line(Picture_ii.FEs(1:inds),Picture_ii.value(1:inds),'Marker',marks_k{ii},'Color',colors{mod(ind,length(Algorithms))+1},'LineStyle','-','LineWidth',1.25,'MarkerSize',7,'Parent',ax2);
                                end
                            end
                        end
                        ind = ind + 1;
                        hold on;
                    end
                    
                    if LM==1 
                        
                        ax1 = gca;
                        ylimits = get(ax1,'YLim');
                        set(ax1,'YLim',1.3*ylimits,'LineWidth',2,'FontSize',17,'box','on');
                        
                        if length(Metric)>1
                            lgd = legend(legends{:},'Location','northwest','NumColumns',2);
                            ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
                            'Color','none', 'XColor','k','YColor','k','LineWidth',2,'FontSize',17);
                        else
                            lgd = legend(legends{:},'Location','northeast','NumColumns',2);
                        end
                        xlabel('FEs','Interpreter','latex','Fontname','Times New Roman','FontSize',17);
                        ylabel(Metric{LM},'Interpreter','latex','Fontname','Times New Roman','FontSize',17);
                        title(lgd,Metric{LM});
                        
                    else
                        ylimits = get(ax2,'YLim');
                        set(ax2,'YLim',1.3*ylimits);
                        xlabel('FEs','Interpreter','latex','Fontname','Times New Roman','FontSize',17);
                        ylabel(Metric{LM},'Interpreter','latex','Fontname','Times New Roman','FontSize',17);
                        lgd = legend(legends{:},'Location','northeast','NumColumns',2);
                        title(lgd,Metric{LM})
                    end
                    
                end
                
%                 xlabel('FEs','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
%                 ylabel(Metric,'Interpreter','latex','Fontname','Times New Roman','FontSize',10);
%                 legend(Algorithms{:});
                
				folder = fullfile('Figures','Convergence',['Comparision_on' sprintf('_%s',legends{:})]);
                [~,~]  = mkdir(folder);
                print(gcf,fullfile(folder,sprintf('%s',Problems_name)),'-depsc');
                
                close all;
            end
        end
        
        function [Data,Result] = CalMetrics(Data,Metrics,ReCMetrics)
            if nargin<3
                ReCMetrics = Metrics;
            end
            Result.runtime = Data.runtime;
            Result.upper_FEs = Data.upper_FEs;
            Result.lower_FEs = Data.lower_FEs;
            for l=length(Data.result)
                if Data.M(1)>1
                    if length(Data.result(l).Population)>Data.N(3)
                        Data.result(l).Population = NSGAII_Update(Data.result(l).Population,Data.N(3),'upper');
                        ReCMetrics = Metrics;
                    end
                end
                Population_l = Data.result(l).Population;
                if ~isempty(Population_l)
                    Feasible     = Population_l.bilevel_feasibles;
                    Population_l = Population_l(Feasible);
                end
                if ~isempty(Population_l) && ~isempty(Data.PF) && Data.M(1)>1
                    Population_l(any(Population_l.upper_objs>repmat(max(Data.PF,[],1),length(Population_l),1),2))=[];
                    if ~isempty(Population_l)
                        [~,ID] = min(pdist2(Population_l.upper_objs,Data.PF),[],2);
                        Index = all(Population_l.upper_objs<Data.PF(ID,:),2);
                        Population_l = Population_l(~Index);
                    end
                end
                
                for p=1:length(ReCMetrics)
                    metricName = ReCMetrics{p};
                    if strcmp(metricName,'LGD')
                        try
                            llFeasible     = Data.result(l).Population.lower_feasibles;
                            Data.result(l).(metricName) = GLOBAL.Metric(str2func(metricName),Data.result(l).Population(llFeasible).upper_best,Data);
                        catch
                            Data.result(l).(metricName) = NaN;
                        end
                    elseif strcmp(metricName,'HV')
                        try
                            Data.result(l).(metricName) = GLOBAL.Metric(str2func(metricName),Data.result(l).Population,Data);
                        catch
                            Data.result(l).(metricName) = NaN;
                        end
                    else
                        try
                            Data.result(l).(metricName) = GLOBAL.Metric(str2func(metricName),Population_l.upper_best,Data);
                        catch
                            Data.result(l).(metricName) = NaN;
                        end
                    end
                end
            end
            for p=1:length(Metrics)
                metricName = Metrics{p};
                Result.(metricName) = Data.result(end).(metricName);
            end
            folder = fullfile('Data',func2str(Data.algorithm));
            save(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',func2str(Data.algorithm),class(Data.problem),Data.D(1),Data.D(2),Data.run)),'Data','Result');
        end
        
        function SaveFig(obj,type,~)
            Pop = [];
            PF = obj(1).PF;
            
            if length(obj)>1
                for i = 1:length(obj)
                    Pop = cat(2,Pop,obj(i).result(end).Population);
                end
                
                PF = obj(1).PF;
                
                [~,ID] = min(pdist2(Pop.upper_objs,obj(end).PF),[],2);
                Index = all(Pop.upper_objs<obj(end).PF(ID,:),2);
                Pop = Pop(~Index);
                
                Pop = Pop.upper_best;
            else
                Pop = obj.result(end).Population;
            end
            if isempty(Pop)
                return
            end
            if nargin > 1
                if nargin ==3
                    Fig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_%s',func2str(obj(1).algorithm),class(obj(1).problem),obj(1).M(1),obj(1).M(2),obj(1).D(1),obj(1).D(2),type),'NumberTitle','off','Position',[520,326,880,652],'Visible','on');
                else
                    Fig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_%s',func2str(obj(1).algorithm),class(obj(1).problem),obj(1).M(1),obj(1).M(2),obj(1).D(1),obj(1).D(2),type),'NumberTitle','off','Position',[0,0,1200,500],'Visible','off');
                end
            else
                Fig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_R%d',func2str(obj(1).algorithm),class(obj(1).problem),obj(1).M(1),obj(1).M(2),obj(1).D(1),obj(1).D(2),obj(1).run),'NumberTitle','off','Position',[0,0,1200,500],'Visible','off');
            end
            
            movegui(Fig,'north');
            
            if nargin < 3
                if size(PF,2)>1
                    subplot1 = subplot(1,2,1,'Parent',Fig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                    set(gcf,'CurrentAxes',subplot1)
                    if ~isempty(PF)
                        if size(PF,1)==1
                            plot(PF(:,1),PF(:,2),'MarkerSize',6,'Marker','d','MarkerSize',10,'LineStyle','none','color',[0,0,0]);
                        else
                            plot(PF(:,1),PF(:,2),'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                        
                        end
                        set(subplot1,'NextPlot','add')
                    end
                    %             Pop = obj.result(end).Population;
                    
                    
                    Pop_Objs = Pop.upper_objs;
                    statuses = Pop.upper_CalObjs & Pop.lower_CalObjs;
                    Pop_Objs_temp = Pop_Objs(statuses,:);
                    feasibles = Pop(statuses).upper_feasibles & Pop(statuses).lower_feasibles;
                    if any(feasibles)
                        plot(Pop_Objs_temp(feasibles,1),Pop_Objs_temp(feasibles,2),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                    end
                    
                    infeasibles = (Pop.upper_CalObjs & ~Pop.upper_feasibles)|(Pop.lower_CalObjs & ~Pop.lower_feasibles);
                    if any(infeasibles)
                        plot(Pop_Objs(infeasibles,1),Pop_Objs(infeasibles,2),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    end
                    
                    lower_infeasibles = (Pop.upper_CalObjs & Pop.upper_feasibles) & ~Pop.lower_CalObjs;
                    if any(lower_infeasibles)
                        plot(Pop_Objs(lower_infeasibles,1),Pop_Objs(lower_infeasibles,2),'MarkerSize',6,'Marker','o','LineStyle','none','color','g');
                    end
                else
                    subplot1 = subplot(1,2,1,'Parent',Fig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                    set(gcf,'CurrentAxes',subplot1)
                    if ~isempty(PF)
                        plot(1:length(Pop),PF(:,1),'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                        set(subplot1,'NextPlot','add')
                    end
                    %             Pop = obj.result(end).Population;
                    index = 1:length(Pop);
                    Pop_Objs = Pop.upper_objs;
                    statuses = Pop.upper_CalObjs & Pop.lower_CalObjs;
                    index_statuses = index(statuses);
                    Pop_Objs_temp = Pop_Objs(statuses,:);
                    feasibles = Pop(statuses).upper_feasibles & Pop(statuses).lower_feasibles;
                    if any(feasibles)
                        plot(index_statuses(feasibles,1),Pop_Objs_temp(feasibles,1),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                    end
                    
                    infeasibles = (Pop.upper_CalObjs & ~Pop.upper_feasibles)|(Pop.lower_CalObjs & ~Pop.lower_feasibles);
                    if any(infeasibles)
                        plot(statuses(infeasibles),Pop_Objs(infeasibles,1),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    end
                    
                    lower_infeasibles = (Pop.upper_CalObjs & Pop.upper_feasibles) & ~Pop.lower_CalObjs;
                    if any(lower_infeasibles)
                        plot(index_statuses(lower_infeasibles),Pop_Objs(lower_infeasibles,1),'MarkerSize',6,'Marker','o','LineStyle','none','color','g');
                    end
                end
                %             set(subplot1,'NextPlot','add')
                %             plot(PF(:,1),PF(:,2),'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                
                xlabel('$F_1$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                ylabel('$F_2$','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                title('Upper Objective Space','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                
                subplot2 = subplot(1,2,2,'Parent',Fig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                plot(transpose(Pop.decs),'Parent',subplot2);
                xlim([1 sum(obj(1).D)]);
                xlabel('Dim','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                ylabel('Value','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                title('Decision Space','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                
            else
                
                if ~isempty(PF)
                    if size(PF,1)==1
                        plot(PF(:,1),PF(:,2),'MarkerSize',6,'Marker','p','MarkerSize',10,'LineStyle','none','color',[0,0,0]);
                    else
                        plot(PF(:,1),PF(:,2),'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
                        
                    end
                    hold on;
                end
                %             Pop = obj.result(end).Population;
                
                if ~isempty(Pop)
                    Feasible     = Pop.bilevel_feasibles;
                    FPop = Pop(Feasible).upper_best;
                end
                if ~isempty(FPop) && ~isempty(PF)
                    [~,ID] = min(pdist2(FPop.upper_objs,PF),[],2);
                    Index = all(FPop.upper_objs<PF(ID,:),2);
                    FPop = FPop(~Index);
                end
                
                if ~isempty(FPop)
                    plot(FPop.upper_objs(1),FPop.upper_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                    legend('A published bilevel optimal solution','The obtained bilevel optimal solutions','FontSize',12);
                end
                
                Pop_Objs = Pop.upper_objs;
%                 statuses = Pop.upper_CalObjs & Pop.lower_CalObjs;
%                 Pop_Objs_temp = Pop_Objs(statuses,:);
%                 feasibles = Pop(statuses).upper_feasibles & Pop(statuses).lower_feasibles;
%                 if any(feasibles)
%                     plot(Pop_Objs_temp(feasibles,1),Pop_Objs_temp(feasibles,2),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
%                 end
                
                infeasibles = (Pop.upper_CalObjs & ~Pop.upper_feasibles)|(Pop.lower_CalObjs & ~Pop.lower_feasibles);
                if any(infeasibles)
                    plot(Pop_Objs(infeasibles,1),Pop_Objs(infeasibles,2),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    if ~isempty(FPop)
                        legend('A published bilevel optimal solution','The obtained bilevel optimal solutions','The obtained infeasible solutions','FontSize',12);
                    else
                        legend('A published bilevel optimal solution','The obtained infeasible solutions','FontSize',12);
                    end
                end
                
%                 lower_infeasibles = (Pop.upper_CalObjs & Pop.upper_feasibles) & ~Pop.lower_feasibles;
%                 if any(lower_infeasibles)
%                     plot(Pop_Objs(lower_infeasibles,1),Pop_Objs(lower_infeasibles,2),'MarkerSize',6,'Marker','o','LineStyle','none','color','g');
%                 end
                
                
            end
            folder = fullfile('Figures',func2str(obj(1).algorithm));
            [~,~]  = mkdir(folder);
            if nargin > 1
%                 saveas(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_%s.fig',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),type)));
                print(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_%s',func2str(obj(1).algorithm),class(obj(1).problem),obj(1).D(1),obj(1).D(2),type)),'-depsc');
%                 saveas(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_%s.svg',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),type)));
            else
                print(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_R%d',func2str(obj(1).algorithm),class(obj(1).problem),obj(1).D(1),obj(1).D(2),obj(1).run)),'-depsc');
%                 saveas(gcf,fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.fig',func2str(obj.algorithm),class(obj.problem),obj.D(1),obj.D(2),obj.run)));
            end
            close('all');
        end
    end
end
