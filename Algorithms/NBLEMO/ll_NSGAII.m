function [Data,Archive,N_Data,T_l_max] = ll_NSGAII(Data,Upper_Data_Pop,Archive,N_A,tao,Epsilon_l,T_l,localsearch)
    
    if nargin < 8
        localsearch = false;
    end
    
    if nargin < 7
        T_l_max = [];
    end
    
    Global = GLOBAL.GetObj();
    
    N_Data = 0;
    
    for i = 1:length(Data)
        
        upper_dec_i = Data(i).upper_dec;
        
        if nargin < 8
%             N_l = N_l_max;
            T_l = inf;
%         else
%             if ~isempty(Archive.Population)
%                 Archive_upper_decs = Archive.Population.upper_decs;
%                 Delta_U = max(max(pdist2(Archive_upper_decs,Archive_upper_decs),[],'all'),eps);
%                 Delta_u = min(pdist2(upper_dec_i,Archive_upper_decs),[],2);
%                 N_l = min(max(round(N_l_max*Delta_u/Delta_U),4),N_l_max);
%                 T_l = min(T_l_max*Delta_u/Delta_U,T_l_max);
%             else
%                 N_l = N_l_max;
%                 T_l = T_l_max;
%             end
        end
        
        lower_Population_i = Data(i).Population;
        N_l = length(lower_Population_i);
        
        FrontNo = Data(i).FrontNo;
        CrowdDis = Data(i).CrowdDis;
        
%         Global.lower_Output(lower_Population_i);
        
        t_i = 0;
        Objs = {};
        
        while t_i < T_l 
            
            try
                Dis = min(pdist2(Archive.Population.upper_decs,upper_dec_i),[],2);
                MatingPop = Archive.Population(Dis==0);
                NDS_MatingPop = NDSort(MatingPop.lower_objs,MatingPop.lower_cons,1);
                CrowdDis_MatingPop = CrowdingDistance(MatingPop.upper_objs,NDS_MatingPop);
                MatingPool = TournamentSelection(2,N_l,NDS_MatingPop,-CrowdDis_MatingPop);
                lower_decs = GA(MatingPop(MatingPool),'lower',{0.9,15,0.1,20});
                lower_Offspring_i = Global.Evaluate(INDIVIDUAL({upper_dec_i,lower_decs}),'lower');
                
            catch
                MatingPool = TournamentSelection(2,N_l,FrontNo,-CrowdDis);
                lower_decs = GA(lower_Population_i(MatingPool),'lower',{0.9,15,0.1,20});
                lower_Offspring_i = Global.Evaluate(INDIVIDUAL({upper_dec_i,lower_decs}),'lower');
                
            end
            
            [lower_Population_i,~,FrontNo,CrowdDis] = NSGAII_Update([lower_Population_i,lower_Offspring_i],N_l,'lower');
            
            Global.lower_Output(lower_Population_i);
            
            Objs = cat(2,Objs,{lower_Population_i(FrontNo==1).lower_objs});
            
            if length(Objs)==tao
                RefPoint = max(cat(1,Objs{:}),[],1);
                HV = cell2mat(cellfun(@(x) CalHV(x,RefPoint), Objs, 'UniformOutput',false));
                H_l = (max(HV)-min(HV))/(max(HV)+min(HV));
                if (H_l <= Epsilon_l && any(lower_Population_i.lower_feasibles))|| t_i>500
                    break;
                end
                Objs(1) = [];
            end
            
            t_i = t_i + 1;
            
        end
        
        lower_Population_i = Global.Evaluate(lower_Population_i,'upper');
        
        [~,best_i] = lower_Population_i.lower_best;
        NDS_best_i = NDSort([lower_Population_i(best_i).upper_objs;Upper_Data_Pop.Population.upper_objs],...
                            [lower_Population_i(best_i).upper_cons;Upper_Data_Pop.Population.upper_cons],1);
        best_i = best_i(NDS_best_i(1:length(best_i))==1);
        
        %% local search operator
        if ~isempty(best_i)
            
            if localsearch
                
                f_max = max(lower_Population_i.lower_objs,[],1);
                f_min = min(lower_Population_i.lower_objs,[],1);
                
                ls_Pop = [];
                
                for ii=best_i
                    
                    f_z = lower_Population_i(ii).lower_obj;
                    
                    options = optimoptions('fmincon','Display','off','Algorithm','sqp','OptimalityTolerance',1e-2);%,'MaxFunctionEvaluations',1e4);
                    problem.options = options;
                    problem.solver = 'fmincon';
                    problem.objective = @(x)CalObj({upper_dec_i,x},f_z,f_max,f_min);
                    problem.nonlcon = @(x)CalCon({upper_dec_i,x});
                    problem.lb = Global.lower_domain(1,:);
                    problem.ub = Global.lower_domain(2,:);
                    problem.x0 = lower_Population_i(ii).lower_dec;
                    
                    [x,~,~,output,~,~,~] = fmincon(problem);
                    Global.SetFEs([Global.upper_FEs,Global.lower_FEs + output.funcCount])
                    
                    Best_ii = Global.Evaluate(INDIVIDUAL({upper_dec_i,x}));
                    Best_ii.label = true;
                    
                    ls_Pop = cat(2,ls_Pop,Best_ii);
                    
                    if ~isempty(Archive.Population)
                        Archive_Population = [Best_ii,Archive.Population];
                        NDS_Archive_Population = NDSort(Archive_Population.upper_objs,Archive_Population.upper_cons,1);
                        CrowdDis_Archive_Population = CrowdingDistance(Archive_Population.upper_objs,NDS_Archive_Population);
                        
                        Archive.Population = Archive_Population(NDS_Archive_Population==1);
                        Archive.CrowdDis = CrowdDis_Archive_Population(NDS_Archive_Population==1);
                        
                        if length(Archive.Population)>N_A
                            [Archive.Population,~,~,Archive.CrowdDis] = NSGAII_Update(Archive.Population,N_A,'lower');
                        end
                    else
                        Archive.Population = Best_ii;
                        Archive.CrowdDis = inf;
                    end
                end
                
                [lower_Population_i,~,FrontNo,CrowdDis] = NSGAII_Update([lower_Population_i,ls_Pop],N_l,'lower');
                
            else
                if ~isempty(Archive.Population)
                    Archive_Population = [lower_Population_i(best_i),Archive.Population];
                    NDS_Archive_Population = NDSort(Archive_Population.upper_objs,Archive_Population.upper_cons,1);
                    CrowdDis_Archive_Population = CrowdingDistance(Archive_Population.upper_objs,NDS_Archive_Population);
                    
                    Archive.Population = Archive_Population(NDS_Archive_Population==1);
                    Archive.CrowdDis = CrowdDis_Archive_Population(NDS_Archive_Population==1);
                    
                    if length(Archive.Population)>N_A
                        [Archive.Population,~,~,Archive.CrowdDis] = NSGAII_Update(Archive.Population,N_A,'lower');
                    end
                else
                    Archive.Population = lower_Population_i(best_i);
                    Archive.CrowdDis = Inf(1,length(best_i));
                end
            end
            
        end
        
        Global.lower_Output(lower_Population_i);
        
        Data(i).upper_dec = upper_dec_i;
        Data(i).Population = lower_Population_i;
        Data(i).FrontNo = FrontNo;
        Data(i).CrowdDis = CrowdDis;
        
        Upper_Data_Pop = Upper_Data(Data);
        
        if nargin < 7
            T_l_max = cat(2,T_l_max,t_i);
        end
        
        N_Data = N_Data + N_l;
        
    end
    
    if nargin < 7
        T_l_max = mean(T_l_max);
    else
        T_l_max = [];
    end
    
end

function Obj = CalObj(x,f_z,f_max,f_min)
    Global = GLOBAL.GetObj();
    Obj = Global.problem.CalObj(x,'lower');
    Obj = max((Obj-f_z)/(f_max-f_min))+1e-6*sum((Obj-f_z)/(f_max-f_min));
end

function [c,ceq] = CalCon(x)
    Global = GLOBAL.GetObj();
    c = Global.problem.CalCon(x,'lower');
    ceq = [];
end