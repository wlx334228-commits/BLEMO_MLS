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
                    if strcmp(err.identifier,'GLOBAL:Termination')
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
            
            % Detect whether the number of evaluations has exceeded
            if nargin < 3
                if any([obj.upper_FEs,obj.lower_FEs] > obj.maxFEs)
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
                                            objs = obj.problem.CalObj(Population(index),'upper');
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
                                            objs = obj.problem.CalObj(Population(index),'lower');
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
                xlim([1 sum(obj.D)]);
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
                    plot(Pop_temp(feasibles).upper_objs(1),Pop_temp(feasibles).upper_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','r');
                    
                    infeasibles = (Pop.upper_CalObjs & ~Pop.upper_feasibles)|(Pop.lower_CalObjs & ~Pop.lower_feasibles);
                    plot(Pop(infeasibles).upper_objs(1),Pop(infeasibles).upper_objs(2),'MarkerSize',6,'Marker','o','LineStyle','none','color','b');
                    
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
                xlim([1 sum(obj.D)]);
                xlabel('Dim','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                ylabel('Value','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                title('Decision Space','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
            end
        end 
        
    end
    
    methods(Access = private, Static)
        %% Display or save the result after the algorithm is terminated
        function Output(obj,notermination)
            if notermination
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
                    if ~isempty(Population_i) && ~isempty(Data.PF)
                        [~,ID] = min(pdist2(Population_i.upper_objs,Data.PF),[],2);
                        Index = all(Population_i.upper_objs<Data.PF(ID,:),2);
                        Population_i = Population_i(~Index);
                    end
                    
                    if length(obj.CalMetric) > 3
                        for j=4:length(obj.CalMetric)
                            metricName = obj.CalMetric{j};
                            if strcmp(metricName,'LGD')
                                try
                                    Data.result(i).(metricName) = GLOBAL.Metric(str2func(metricName),Data.result(i).Population.upper_best,Data);
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
        
        function Data_dealing(Problems,Algorithms,Run,Metrics,Re,type)
            metric_names = {'runtime','upper_FEs','lower_FEs','FEs'};
            if nargin > 3
                metric_names = cat(2,metric_names,Metrics);
            end
            
            if nargin > 4
                RecalM = Re;
            else
                RecalM = false;
            end
            
            if nargin<5
                type = 'std';
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
                    Population = [];
                    
                    for k=1:Run
                        
                        load(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',Algorithms{i},Problems_name,ProblemsD(1),ProblemsD(2),k)));
                        
                        if size(Data.PF,1)<=1 && Data.M(1)>1
                            RecalM = true;
                            Population = cat(2,Population,Data.result(end).Population);
                        end
                        
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
                    end
                    
                    if RecalM

                        if ~isempty(Population)
                            PF = [Data.problem.PF();Population.upper_best.upper_objs];
                        else
                            try 
                                PF = Data.problem.PF;
                            catch
                                [~,~,PF,~]=getOptimalSolutionSMD(Data.D(1),Data.D(2),func2str(Data.problem));
                            end
                        end

                        for k=1:Run
                            
                            Global(k).PF = PF;
                            
                            [Global(k),Result] = GLOBAL.CalMetrics(Global(k),Metrics);
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
                        end
                        
                        switch type
                            case 'std'
                                Vars = std(result_k(ID));
                            case 'iqr'
                                Vars = iqr(result_k(ID));
                        end
                        
                        if ~isempty(ID)
                            if length(Algorithms)>1 && i>1 && k>4
                                ranksum_result_j = ranksum_result.(Problems{j,1});
                                ranksum_result_jk = ranksum_result_j.(metric_names{k});
                                
                                [~,p] = ranksum(result_k(ID),ranksum_result_jk);
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
                                    mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)%s',mean(result_k(ID)),Vars,'=');
                                    %                                     mean_result.(metric_names{k}) =  sprintf('%0.3d%s',mean(result_k(ID)),'=');
                                end
                            else
                                if k>4
                                    mean_result.(metric_names{k}) =  sprintf('%0.2e(%0.2e)',mean(result_k(ID)),Vars);
                                else
                                    mean_result.(metric_names{k}) = mean(result_k(ID));
                                end
                            end
                        else
                            if length(Algorithms)>1 && i>1 && k>4
                                mean_result.(metric_names{k}) = 'NaN(NaN)+';
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
        
        function Convergence_plot(Problems,Algorithms,Metric,Runs)
            close all;
            for j=1:size(Problems,1)
                
                Problems_name = Problems{j,1};
                ProblemsD = Problems{j,2};
                figure('Name',sprintf('The Convergence Plot on %s',Problems_name),'Position',[520,343,880,635],'Visible','on');
                marks ={'o','d','s','*','p','h'};
                colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#A2142F'};
                legends = {'BLEA-UDP','SMS-EMOA','H-BLEMO','NBLEMO'};
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
                    
                    for k=1:Runs
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
                                    UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),[Metric_k,zeros(1,Y-length(Metric_k))]);
                                    FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),[FE_k,zeros(1,Y-length(FE_k))]);
                                else
                                    UMetric.(Metric{LM})(:,Y+1:length(Metric_k)) = 0;
                                    FE.(Metric{LM})(:,Y+1:length(FE_k)) = 0;
                                    UMetric.(Metric{LM}) = cat(1,UMetric.(Metric{LM}),Metric_k);
                                    FE.(Metric{LM}) = cat(1,FE.(Metric{LM}),FE_k);
                                end
                            end
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
                        
                        II = FEs>0 & ~isnan(UMetric_LM);
                        
                        I1 = any(II==true,2);
                        I2 = any(II==true,1);
                        
                        II = II(I1,I2);
                        sumN = sum(II,1);
                        FEs = FEs(I1,I2);
                        FEs(~II) = 0;
                        FEs = sum(FEs,1)./(sumN+eps);
                        
                        UMetric_LM = UMetric_LM(I1,I2);
                        UMetric_LM(isnan(UMetric_LM))=0;
                        UMetric_LM = sum(UMetric_LM,1)./(sumN+eps);
                        
                        [~,rank] =sort(FEs,'ascend');
                        FEs = FEs(rank);
                        UMetric_LM = UMetric_LM(rank);
                        
                        if strcmp(Metric{LM},'UHV')
                            DV = diff(UMetric_LM);
                            index = DV<0 & abs(DV)>0.01;
                            while any(index)
                                index = find(index)+1;
                                UMetric_LM(index)=[];
                                FEs(index)=[];
                                DV = diff(UMetric_LM);
                                index = DV<0 & abs(DV)>0.01;
                            end
                        else
                            DV = diff(UMetric_LM);
                            index = DV>0 & abs(DV)>0.01;
                            while any(index)
                                index = find(index)+1;
                                UMetric_LM(index)=[];
                                FEs(index)=[];
                                DV = diff(UMetric_LM);
                                index = DV>0 & abs(DV)>0.01;
                            end
                        end
                        
                        Picture_i.(Metric{LM}).FEs = [0,FEs];
                        Picture_i.(Metric{LM}).maxFEs = max(FEs);
                        if strcmp(Metric{LM},'UHV')
                            Picture_i.(Metric{LM}).value = [0,UMetric_LM];
                        else
                            Picture_i.(Metric{LM}).value = [1.1,UMetric_LM];
                        end
                    end
                    
                    Picture = cat(1,Picture,Picture_i);
                    
                end
                
                ind = 1;
                
                for LM = 1:length(Metric)
                    Picture_i = [Picture.(Metric{LM})];
                    switch Problems_name(1:2)
                        case 'TP'
                            maxFEs = 1.2*min([Picture_i.maxFEs]);
                        otherwise
                            maxFEs = min([Picture_i.maxFEs]);
                    end
                    
                    
                    for ii=1:length(Picture_i)
                        Picture_ii = Picture_i(ii);
                        inds = find(Picture_ii.FEs<=maxFEs, 1, 'last' );
                        if inds<length(Picture_ii.FEs)
                            if LM == 1
                                line([Picture_ii.FEs(1:inds),maxFEs],Picture_ii.value(1:inds+1),'Marker',marks{ii},'Color',colors{ii},'LineWidth',1,'MarkerSize',5);
                            else
                                line([Picture_ii.FEs(1:inds),maxFEs],Picture_ii.value(1:inds+1),'Marker',marks{ii},'Color',colors{ii},'LineStyle','-','LineWidth',1,'MarkerSize',5,'Parent',ax2);
                            end
                        else
                            if LM == 1
                                line(Picture_ii.FEs(1:inds),Picture_ii.value(1:inds),'Marker',marks{ii},'Color',colors{ii},'LineWidth',1,'MarkerSize',5);
                            else
                                line(Picture_ii.FEs(1:inds),Picture_ii.value(1:inds),'Marker',marks{ii},'Color',colors{ii},'LineStyle','-','LineWidth',1,'MarkerSize',5,'Parent',ax2);
                            end
                        end
                        
                        ind = ind + 1;
                        hold on;
                    end
                    
                    if LM==1
                        
                        ax1 = gca;
                        ylimits = get(ax1,'YLim');
                        set(ax1,'YLim',1.15*ylimits);
                        xlabel('FEs','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                        ylabel(Metric{LM},'Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                        lgd = legend(legends{:},'Location','northwest','NumColumns',2);
                        title(lgd,Metric{LM});
                        
                        ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
                            'Color','none', 'XColor','k','YColor','k');
                    else
                        ylimits = get(ax2,'YLim');
                        set(ax2,'YLim',1.5*ylimits);
                        xlabel('FEs','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                        ylabel(Metric{LM},'Interpreter','latex','Fontname','Times New Roman','FontSize',10);
                        lgd = legend(legends{:},'Location','northeast','NumColumns',2);
                        title(lgd,Metric{LM})
                    end
                    
                end
                
%                 xlabel('FEs','Interpreter','latex','Fontname','Times New Roman','FontSize',10);
%                 ylabel(Metric,'Interpreter','latex','Fontname','Times New Roman','FontSize',10);
%                 legend(Algorithms{:});
                
                folder = 'Figures\Convergence';
                [~,~]  = mkdir(folder);
                print(gcf,fullfile(folder,sprintf('%s',Problems_name)),'-depsc');
                
                close all;
            end
        end
        
        function [Data,Result] = CalMetrics(Data,Metrics)
            Result.runtime = Data.runtime;
            Result.upper_FEs = Data.upper_FEs;
            Result.lower_FEs = Data.lower_FEs;
            for l=1:length(Data.result)
                Population_l = Data.result(l).Population;
                Feasible     = Population_l.bilevel_feasibles;
                Population_l = Population_l(Feasible);
%                 if ~isempty(Population_l)
%                     [~,ID] = min(pdist2(Population_l.upper_objs,Data.PF),[],2);
%                     Index = all(Population_l.upper_objs<Data.PF(ID,:),2);
%                     Population_l = Population_l(~Index);
%                 end
                for p=1:length(Metrics)
                    metricName = Metrics{p};
                    Data.result(l).(metricName) = GLOBAL.Metric(str2func(metricName),Population_l.upper_best,Data);
                    if l==length(Data.result)
                        Result.(metricName) = Data.result(l).(metricName);
                    end
                end
            end
            folder = fullfile('Data',func2str(Data.algorithm));
            save(fullfile(folder,sprintf('%s_%s_D%d_%d_R%d.mat',func2str(Data.algorithm),class(Data.problem),Data.D(1),Data.D(2),Data.run)),'Data','Result');
        end
        
        function SaveFig(obj,type)
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
            
            if nargin > 1
                Fig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_%s',func2str(obj(1).algorithm),class(obj(1).problem),obj(1).M(1),obj(1).M(2),obj(1).D(1),obj(1).D(2),type),'NumberTitle','off','Position',[0,0,1200,500],'Visible','off');    
            else
                Fig = figure('Name',sprintf('%s_%s_M%d_%d_D%d_%d_R%d',func2str(obj(1).algorithm),class(obj(1).problem),obj(1).M(1),obj(1).M(2),obj(1).D(1),obj(1).D(2),obj(1).run),'NumberTitle','off','Position',[0,0,1200,500],'Visible','off');
            end
            
            movegui(Fig,'north');
            if size(PF,2)>1
                subplot1 = subplot(1,2,1,'Parent',Fig,'Fontname','Times New Roman','FontSize',10,'box','on','NextPlot','replacechildren');
                set(gcf,'CurrentAxes',subplot1)
                if ~isempty(PF)
                    plot(PF(:,1),PF(:,2),'MarkerSize',6,'Marker','.','LineStyle','none','color',[0,0,0]);
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
                
                check_infeasibles = (Pop.upper_CalObjs & Pop.upper_feasibles) & ~Pop.lower_CalObjs;
                if any(check_infeasibles)
                    plot(Pop_Objs(check_infeasibles,1),Pop_Objs(check_infeasibles,2),'MarkerSize',6,'Marker','o','LineStyle','none','color','g');
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
               
               check_infeasibles = (Pop.upper_CalObjs & Pop.upper_feasibles) & ~Pop.lower_CalObjs;
               if any(check_infeasibles)
                   plot(index_statuses(check_infeasibles),Pop_Objs(check_infeasibles,1),'MarkerSize',6,'Marker','o','LineStyle','none','color','g');
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