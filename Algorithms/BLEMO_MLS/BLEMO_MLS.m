function BLEMO_MLS(Global)

    if Global.N(1) <=10
        [K,T_u,T_l,alpha,belta,theta,LMaxG] = Global.ParameterSet(5,5,5,10,0.1,1e-3,500);
    else
        [K,T_u,T_l,alpha,belta,theta,LMaxG] = Global.ParameterSet(10,10,5,10,0.1,1e-3,500);
    end
    
    %Q函数更新参数
    pha0 = 0.9;
    lmd = 0.8;
    
    %初始化Q表
    Actions = 1:3;
    Status = 1:4;
    Qtable = zeros(length(Status),length(Actions));
    hvWindow = 10;
    hvTolerance = 1e-4;
    minUpperGen = 100;
    
    W = UniformPoint(K,Global.M(1));
    N = Global.N(3);
    upper_N = Global.N(1);
    History_data = [];
    Elite_Population = [];
    
    [redundant_dec_uobj, redundant_dec_lobj] = group(Global); %决策变量分组
    
    Population_A = Global.Evaluate(Global.Initialization);
    [Population_B,History_data] = lower_level_optimize(Population_A,History_data,redundant_dec_lobj,T_l,alpha,belta,LMaxG);
    NDrank = [];
    for i=1:length(Population_A)
        Population_i = [Population_A(i),Population_B(i)];
        NDrank = cat(1,NDrank,NDSort(Population_i.upper_objs,Population_i.upper_cons,inf));
    end
    
    if sum(NDrank(:,1)>=NDrank(:,2))/upper_N>=0.9 % 判断问题类型，第一类问题进行忽略下层制约的上层探索
        Population_A = tentative_explore(Population_B,redundant_dec_uobj,alpha,theta);
        [Population_A,History_data] = lower_level_optimize(Population_A,History_data,redundant_dec_lobj,T_l,alpha,belta,LMaxG);
        Population = NSGAII_Update([Population_A,Population_B],upper_N,'upper');
    else
        Population = Population_B;
    end
    
    clear Population_A Population_B Population_i NDrank
    
    [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,redundant_dec_lobj,W,T_l,alpha,belta,LMaxG,1);
    upperHVArchive = UpdateTerminationArchive({},UpperTerminationObjs(Elite_Population),hvWindow);
    
    Global.upper_Output(Population,2);
    
    Last_Ac = 1;
    current_Ac = randsample(Actions,1);
    
    Erate = 1;
    Rrate = 1;
    current_St = 1;
    
    while Global.NotTermination(Elite_Population,UpperTerminationContinue(Global,upperHVArchive,hvTolerance,minUpperGen))
        
        ID = find(Population.labels);
        [Population(ID),History_data] = lower_level_optimize(Population(ID),History_data,redundant_dec_lobj,T_l,alpha,belta,LMaxG);
        [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,redundant_dec_lobj,W,T_l,alpha,belta,LMaxG,1);
        
        Offspring = Generate(Population,History_data,redundant_dec_uobj,W,T_u);
        Population = cat(2,Population,Offspring);
        
        [~,History_Population] = Update_Data(History_data,upper_N);
        
        % 执行选定下层搜索模式
        switch current_Ac
            case 1 %对所有未经下层搜索的个体执行LLS
                ID = ~Population.labels;
                [Population(ID),History_data] = lower_level_optimize(Population(ID),History_data,redundant_dec_lobj,T_l,alpha,belta,LMaxG);
                
                %rate of effective LLS
                if ~isempty(ID)
                    I_epsilon = Cal_I_epsilon(Population(ID).upper_objs,History_Population.upper_objs);
                    New_Erate = sum(I_epsilon<=1 & Population(ID).upper_feasibles)/length(I_epsilon);
                else
                    New_Erate = Erate;
                end
                
                Population = Upper_Select(Population,upper_N,W,'upper');
                
            case 2 % 仅对有前景的个体执行LLS
                ID = Candidate_choose(Population,upper_N,W);
                [Population(ID),History_data] = lower_level_optimize(Population(ID),History_data,redundant_dec_lobj,T_l,alpha,belta,LMaxG);
                
                %rate of effective LLS
                if ~isempty(ID)
                    I_epsilon = Cal_I_epsilon(Population(ID).upper_objs,History_Population.upper_objs);
                    New_Erate = sum(I_epsilon<=1 & Population(ID).upper_feasibles)/length(I_epsilon);
                else
                    New_Erate = Erate;
                end
                
                Population = Upper_Select(Population,upper_N,W,'upper');
                
            case 3 %不执行LLS

                Population = Upper_Select(Population,upper_N,W,'upper');
                
                New_Erate = Erate;
        end
        
        % 更新档案
        History_data = Update_Data(History_data,10*N);
        [Elite_Population,History_data] = Update_Elites(Elite_Population,History_data,redundant_dec_lobj,W,T_l,alpha,belta,LMaxG,2);

        %% Update status
        NDPop = Population.upper_best;
        New_Rrate = sum(NDPop.labels)/length(NDPop);
        
        St_E = New_Erate >= Erate;
        St_R = New_Rrate >= Rrate;
        
        switch current_Ac
            case {1 2}
                    if St_E
                        New_St = 1;
                    else
                        New_St = 2;
                    end
            case 3
                switch current_St
                    case {1 2 3}
                        if St_R && New_Rrate>0
                            New_St = 3;
                        else
                            New_St = 4;
                        end
                    case 4
                        if St_R
                            New_St = 3;
                        else
                            New_St = 4;
                        end
                end
        end
        
        %% Update Qtable
        progress = min(TotalFE(Global)/max(TotalMaxFE(Global),eps),1);
        pha = 1-(pha0*progress);
        
        if current_Ac~=3 || (current_Ac==3 && current_St==4)
            if New_St<3
                if New_St==1
                    Qtable(current_St,current_Ac) = Qtable(current_St,current_Ac)+pha*(max(0.01,New_Erate)+ lmd*max(Qtable(New_St,:))-Qtable(current_St,current_Ac));
                else
                    Qtable(current_St,current_Ac) = Qtable(current_St,current_Ac)+pha*(-max(0.01,abs(New_Erate-Erate))+ lmd*max(Qtable(New_St,:))-Qtable(current_St,current_Ac));
                end
            else
                if New_St==3
                    Qtable(current_St,current_Ac) = Qtable(current_St,current_Ac)+pha*(max(0.01,New_Rrate)+ lmd*max(Qtable(New_St,:))-Qtable(current_St,current_Ac));
                else
                    Qtable(current_St,current_Ac) = Qtable(current_St,current_Ac)+pha*(-max(0.01,abs(New_Rrate-Rrate))+ lmd*max(Qtable(New_St,:))-Qtable(current_St,current_Ac));
                end
            end
        end
        
        if Last_Ac ==3 && Last_St~=4
            switch current_Ac
                case {1 2}
                    if New_St==1
                        Qtable(Last_St,Last_Ac)  = Qtable(Last_St,Last_Ac) +pha*(max(0.01,New_Erate)+ lmd*max(Qtable(New_St,:))-Qtable(Last_St,Last_Ac));
                    else
                        Qtable(Last_St,Last_Ac)  = Qtable(Last_St,Last_Ac) +pha*(-max(0.01,abs(New_Erate-Erate))+ lmd*max(Qtable(New_St,:))-Qtable(Last_St,Last_Ac));
                    end
                case 3
                    Qtable(Last_St,Last_Ac)  = Qtable(Last_St,Last_Ac) +pha*(-max(0.01,abs(New_Rrate-Rrate))+ lmd*max(Qtable(New_St,:))-Qtable(Last_St,Last_Ac));
            end
        end
        
        Erate = New_Erate;
        Rrate = New_Rrate;
        Last_St = current_St;
        current_St = New_St;
        
        %% Update action
        rnd = 0.5*(1+progress);
        Last_Ac = current_Ac;
        if length(unique(Qtable(current_St,:)))>1 && rand < rnd && all(Qtable(current_St,:)~=0)
            Actions_temp = find(Qtable(current_St,:)==max(Qtable(current_St,:)));
            current_Ac = Actions_temp(randi(length(Actions_temp),1));
        else
            if any(Qtable(current_St,:)==0)
                Actions_temp = Actions(Qtable(current_St,:)==0);
            else
                Actions_temp = Actions;
            end
            if ~isempty(Actions_temp)
                current_Ac = Actions_temp(randi(length(Actions_temp),1));
            else
                current_Ac = randsample(Actions,1);
            end
        end
        
        Global.upper_Output(Population,2); 
        upperHVArchive = UpdateTerminationArchive(upperHVArchive,UpperTerminationObjs(Elite_Population),hvWindow);

        if ~UpperTerminationContinue(Global,upperHVArchive,hvTolerance,minUpperGen)
            Global.NotTermination(Elite_Population,false);
        end
    end
end

function keepGoing = UpperTerminationContinue(Global,hvArchive,hvTolerance,minUpperGen)
    keepGoing = TotalFE(Global) < TotalMaxFE(Global) && ...
        (Global.gen < minUpperGen || ~HVArchiveConverged(hvArchive,hvTolerance));
end

function totalFE = TotalFE(Global)
    totalFE = Global.upper_FEs + Global.lower_FEs;
end

function maxFE = TotalMaxFE(Global)
    if isscalar(Global.maxFEs)
        maxFE = Global.maxFEs;
    else
        maxFE = sum(Global.maxFEs);
    end
end

function hvArchive = UpdateTerminationArchive(hvArchive,PopObj,windowSize)
    if ~isempty(PopObj)
        hvArchive = cat(2,hvArchive,{PopObj});
        if length(hvArchive) > windowSize
            hvArchive(1:length(hvArchive)-windowSize) = [];
        end
    end
end

function PopObj = UpperTerminationObjs(Population)
    PopObj = [];
    if isempty(Population)
        return;
    end
    Feasible = Population.upper_feasibles & Population.lower_feasibles;
    Population = Population(Feasible);
    if isempty(Population)
        return;
    end
    FrontNo = NDSort(Population.upper_objs,1);
    PopObj = Population(FrontNo==1).upper_objs;
end

function converged = HVArchiveConverged(hvArchive,hvTolerance)
    converged = false;
    if length(hvArchive) < 10
        return;
    end
    PopObjAll = cat(1,hvArchive{:});
    if isempty(PopObjAll)
        return;
    end
    RefPoint = max(PopObjAll,[],1) + 0.1;
    HV = cellfun(@(PopObj)TerminationHV(PopObj,RefPoint),hvArchive);
    relHV = (max(HV)-min(HV))/(max(HV)+min(HV)+eps);
    converged = relHV < hvTolerance;
end
