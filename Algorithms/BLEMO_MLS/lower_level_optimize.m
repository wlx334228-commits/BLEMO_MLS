function [Population,History_data] = lower_level_optimize(Population,History_data,redundant_dec,T,alpha,~,Maxgen)
    if isempty(Population)
        return;
    end

    Global = GLOBAL.GetObj();
    
    redundant_dec_l = redundant_dec(redundant_dec>Global.D(1))-Global.D(1);
    
    Population = Global.Evaluate(Population);
    
    N = length(Population);
    
    upper_decs = Population.upper_decs;
    [unique_upper_decs,IA] = uniquetol(upper_decs,0,'ByRows', true,'OutputAllIndices', true);
    if ~isempty(History_data)
        [Dis,Ind] = min(pdist2(unique_upper_decs,History_data.upper_decs),[],2);
    else
        Dis = Inf(1,length(IA));
    end
    
    Archive = cell(length(IA),1);
    
    lower_RV = Population.lower_objs; %指定下层搜索方向
    
    ArchiveID = 1:N;
    
   
    for i=1:length(IA)
        
        ArchiveID(IA{i}) = i;
        
        Population_i = Population(IA{i});
        
        if length(IA{i}) > 1
            [~,I] = unique(Population_i.lower_decs,'row');
            Population_i = Population_i(I);
        end
        
        index = 1:N;
        if isempty(History_data)||Dis(i)>0
            
            IA_i = transpose(IA{i});
            
            index(IA_i) = [];
            if ~isempty(index)
                if length(index)>=Global.N(2)
                    lower_decs = Population(index(randperm(length(index),Global.N(2)))).lower_decs;
                    Archive_i = cat(2,Population_i,INDIVIDUAL({unique_upper_decs(i,:),lower_decs}));
                else
                    lower_decs = Population(index).lower_decs;
                    Archive_i = cat(2,Population_i,INDIVIDUAL({unique_upper_decs(i,:),lower_decs}));
                end
            else
                Archive_i = Population_i;
            end
            
            while length(Archive_i)<Global.N(2)
                lower_decs = mutation(Population_i.lower_decs,1,'lower',{1,20});
                Archive_i = cat(2,Archive_i,INDIVIDUAL({unique_upper_decs(i,:),lower_decs}));
                [~,I] = unique(Archive_i.lower_decs,'row');
                Archive_i = Archive_i(I);
            end
            
            Archive_i = Global.Evaluate(Archive_i,'lower');
            
        else
            Archive_i = History_data(Ind(i)).Archive;
            D = min(pdist2(Population_i.lower_decs,Archive_i.lower_decs),[],2);
            Archive_i = cat(2,Archive_i,Population_i(D~=0));
        end
        
        IA_i = transpose(IA{i});
        
        %% 基于内点法的下层LLS
        for j=IA_i
            
            Population_j = Population(j);
            Zmin = min(Archive_i.lower_objs,[],1);
            Zmin = min(Zmin,0);
            lower_RV_i = lower_RV(i,:);
            Z = max(lower_RV_i-Zmin,eps);
            Z = Z./sum(Z,2);
            
            options = optimoptions('fmincon','Display','off','Algorithm','interior-point','ConstraintTolerance',0,'StepTolerance',1e-2);
            problem.options = options;
            problem.solver = 'fmincon';
            problem.objective = @(x)CalObj({Population_j,redundant_dec_l,x},Zmin,Z);
            problem.nonlcon = @(x)CalCon({Population_j,redundant_dec_l,x});
            problem.lb = Global.lower_domain(1,:);
            problem.ub = Global.lower_domain(2,:);
            problem.x0 = Population_j.lower_dec;
            
            [lower_dec,~,~,output,~,~,~] = fmincon(problem);
            lower_dec(:,redundant_dec_l) = Population_j.lower_dec(redundant_dec_l);
            
            Global.SetFEs([Global.upper_FEs,Global.lower_FEs + output.funcCount])
            
            if min(pdist2(lower_dec,Archive_i.lower_decs),[],2)>1e-6
                Offspring_i = Global.Evaluate(INDIVIDUAL({Population_j.upper_dec,lower_dec}),'lower');
                Archive_i = cat(2,Archive_i,Offspring_i);
                Candidate = [Population_j,Offspring_i];
                
                if Offspring_i.lower_CV == Population_j.lower_CV
                    Population_j = Tchebycheff_Select(Candidate,Archive_i,lower_RV_i,'lower');
                elseif Offspring_i.lower_CV < Population_j.lower_CV
                    Population_j = Offspring_i;
                end
                
                Population(j) = Population_j;
            end
        end
        
        Archive_i = Global.Evaluate(Archive_i,'lower');
        Archive_i = Update_Archive(Archive_i,Global.N(2));
        Archive(i) = {Archive_i};
        
    end
    
    if ~isempty(History_data)
        History_data(Ind(Dis==0))=[];
    end
    
%     Global.lower_Output(Population,1);
    
    t = 0;
    record_index = [];
    indicator = [];
    status = [];
    hvWindow = 10;
    hvTolerance = 1e-3;
    lowerHVArchives = cell(length(Archive),1);
    
    rate = ones(N,1);
    
    %% 局部LLS后的下层进化优化
    while t < Maxgen
        indicator_t = zeros(N,1);
        status_t = zeros(N,1);
        for i = 1:N
            Population_i = Population(i);
            lower_RV_i = lower_RV(i,:);
            if ~ismember(i,record_index) || ...
                    ~Population_i.lower_feasible ||...
                    (ismember(i,record_index) && Population_i.lower_feasible && rand<1-rate(i))
                
                Archive_i = Archive{ArchiveID(i)};
                
                if length(IA) >= T
                    if rand < 0.5 || N < T
                        P = Mating(Population_i,Archive_i,T);
                        lower_dec = DE_current_1_bin(Population_i.lower_dec,P(1).lower_dec,P(2).lower_dec,...
                            'lower',{0.5,0.2,0,20});
                    else
                        P = Mating(Population_i,Population,T);
                        lower_dec = DE_current_1_bin(Population_i.lower_dec,P(1).lower_dec,P(2).lower_dec,...
                            'lower',{0.5,0.2,0,20});
                    end
                else
                    lower_dec = mutation(Population_i.lower_dec,1,'lower',{1,20});
                end
                
                lower_dec(:,redundant_dec_l) = Population_i.lower_dec(:,redundant_dec_l);
                
                Dis = min(pdist2(lower_dec,Archive_i.lower_decs),[],2);
                
                if Dis>1e-6
                    Offspring_i = Global.Evaluate(INDIVIDUAL({Population_i.upper_dec,lower_dec}),'lower');
                      
                    Archive_i = cat(2,Archive_i,Offspring_i);
                    Candidate = [Population_i,Offspring_i];
                    
                    if Offspring_i.lower_CV == Population_i.lower_CV
                        [Population_i,pos] = Tchebycheff_Select(Candidate,Archive_i,lower_RV_i,'lower');
                        if pos==2 && Offspring_i.lower_CV ==0
                            indicator_t(i) = 1;
                        end
                        
                    elseif Offspring_i.lower_CV < Population_i.lower_CV
                        Population_i = Offspring_i;
                        indicator_t(i) = 1;
                    end
                    
                    Archive_i = Update_Archive(Archive_i,Global.N(2));
                    
                    Archive(ArchiveID(i)) = {Archive_i};
                    Population(i) = Population_i;
                    status_t(i) = 1;
                    
                end
            end
        end
        
        indicator = cat(2,indicator,indicator_t);
        status = cat(2,status,status_t);
        
        rate = sum(indicator==1,2)./sum(status,2);
        lowerHVArchives = UpdateLowerTerminationArchives(lowerHVArchives,Archive,hvWindow);
        
        %% 下层LLS终止判定
        if  t>=alpha
            record_index = find(all(indicator(:,t-alpha+1:end)==0,2));
        end
        if LowerTerminationConverged(lowerHVArchives,hvTolerance)
            break;
        end
        
%         Global.lower_Output(Population,3);
        
        t=t+1;
        
    end
    
    Current_data = Data_dealing(unique_upper_decs,Archive);
    Population = Update(Population,Current_data,IA,lower_RV);
    
    Global.lower_Output(Population,3);
    
    History_data = cat(2,History_data,Current_data);
end

function P = Mating(Individual,Data,T)
    lower_obj = Individual.lower_obj;
    Findex = find(Data.lower_feasibles);
    if length(Findex)>T
        B = pdist2(lower_obj,Data(Findex).lower_objs);
        B(B==0) = inf;
        [~,B] = sort(B,2);
        B = B(:,1:T);
        if rand < 0.5
            P = B(randperm(T,2));
        else
            P = randperm(length(Findex),2);
        end
        P = Findex(P);
    else
        FrontNo = NDSort(Data.lower_objs,Data.lower_CVs,inf);
        P=ones(1,2);
        while P(1)==P(2)
            P = TournamentSelection(2,2,FrontNo);
        end
    end
    P = Data(P);
end

function [Population,index] = Tchebycheff_Select(Candidate,Archive,RV,type)
    
    switch type
        case 'upper'
            Zmin = min([RV;Archive.upper_objs],[],1);
        case 'lower'
            Zmin = min([RV;Archive.lower_objs],[],1);
    end
    Zmin = min(Zmin,0);
    RV = max(RV-Zmin,eps);
%     RV = RV-lower_Zmin;
    RV = RV./sum(RV,2);
    index_c = false(1,length(Candidate));
    index = [];
    for i=1:size(RV,1)
        inds = find(~index_c);
        if isempty(inds)
            index_c = false(1,length(Candidate));
            inds = find(~index_c);
        end
        switch type
            case 'upper'
                [~,pos] = min(max(abs(Candidate(inds).upper_objs-Zmin)./max(RV(i,:),eps),[],2));
            case 'lower'
                [~,pos] = min(max(abs(Candidate(inds).lower_objs-Zmin)./max(RV(i,:),eps),[],2));
        end
        index_c(inds(pos)) = true;
        index = cat(2,index,inds(pos));
    end
    Population = Candidate(index);
end

function Archive = Update_Archive(Archive,N)
    if length(Archive)> N
        LFS = find(Archive.lower_feasibles);
        if length(LFS)> N
            PopObj = Archive(LFS).lower_objs;
            [FrontNo,MaxFNo] = NDSort(PopObj,N);
            CrowdDis = CrowdingDistance(PopObj,FrontNo);
            if MaxFNo > 1
                IndexP = FrontNo < MaxFNo;
                Last     = find(FrontNo==MaxFNo);
                [~,Rank] = sort(CrowdDis(Last),'descend');
                IndexP(Last(Rank(1:N-sum(IndexP)))) = true;
            else
                IndexP = Archive(LFS).upper_CalObjs & FrontNo==1;
                if sum(IndexP)<=N
                    Last     = find(FrontNo==MaxFNo & ~IndexP);
                    [~,Rank] = sort(CrowdDis(Last),'descend');
                    IndexP(Last(Rank(1:N-sum(IndexP)))) = true;
                else 
                    Inds = find(IndexP);
                    PopObj = Archive(LFS(Inds)).upper_objs;
                    [FrontNo,MaxFNo] = NDSort(PopObj,Archive(LFS(Inds)).upper_cons,N);
                    IndexP = FrontNo < MaxFNo;
                    Last     = find(FrontNo==MaxFNo);
                    [~,Rank] = sort(CrowdDis(Last),'descend');
                    IndexP(Last(Rank(1:N-sum(IndexP)))) = true;
                    IndexP = Inds(IndexP);
                end
            end
            Archive = Archive(LFS(IndexP));
        else
            [~,rank] = sort(Archive.lower_CVs,'ascend');
            Archive = Archive(rank(1:N));
        end
    end
end

function Data = Data_dealing(upper_decs,Archive)
    Global = GLOBAL.GetObj();
    S = length(Archive);
    Best_Archive = cell(S,1);
    for i = 1:S
        Data_i = Archive{i};
        [~,index] = Data_i.lower_best;
        
        Data_i(index) = Global.Evaluate(Data_i(index),'upper');
        inds = find(Data_i(index).upper_feasibles);
        if ~isempty(inds)
            [~,best] = Data_i(index(inds)).upper_best;
            Best_Archive_i = Data_i(index(inds(best)));
        else
            [~,best] = min(Data_i(index).upper_CVs);
            Best_Archive_i = Data_i(index(best));
        end
        
        for j=1:length(Best_Archive_i)
            Best_Archive_i(j).label = true;
        end
        Archive{i} = Data_i;
        Best_Archive{i} = Best_Archive_i;
    end
    
    Data = DATA(upper_decs,Archive,Best_Archive);
end

function Population = Update(Population,Data,IA,upper_RV)
    S = length(Data);
    for i = 1:S
        Best_Archive_i = Data(i).Best_Archive;
        
        Population(IA{i}) = Tchebycheff_Select(Best_Archive_i,Data(i).Archive,upper_RV(IA{i},:),'lower');
        
    end
    
end

function lowerHVArchives = UpdateLowerTerminationArchives(lowerHVArchives,Archive,windowSize)
    for i = 1:length(Archive)
        PopObj = LowerTerminationObjs(Archive{i});
        if ~isempty(PopObj)
            lowerHVArchives{i} = cat(2,lowerHVArchives{i},{PopObj});
            if length(lowerHVArchives{i}) > windowSize
                lowerHVArchives{i}(1:length(lowerHVArchives{i})-windowSize) = [];
            end
        end
    end
end

function PopObj = LowerTerminationObjs(Archive)
    PopObj = [];
    if isempty(Archive)
        return;
    end
    Archive = Archive(Archive.lower_feasibles);
    if isempty(Archive)
        return;
    end
    FrontNo = NDSort(Archive.lower_objs,1);
    PopObj = Archive(FrontNo==1).lower_objs;
end

function converged = LowerTerminationConverged(lowerHVArchives,hvTolerance)
    converged = false;
    if isempty(lowerHVArchives)
        return;
    end
    for i = 1:length(lowerHVArchives)
        hvArchive = lowerHVArchives{i};
        if length(hvArchive) < 10
            continue;
        end
        PopObjAll = cat(1,hvArchive{:});
        if isempty(PopObjAll)
            continue;
        end
        RefPoint = max(PopObjAll,[],1) + 0.1;
        HV = cellfun(@(PopObj)TerminationHV(PopObj,RefPoint),hvArchive);
        relHV = (max(HV)-min(HV))/(max(HV)+min(HV)+eps);
        if relHV < hvTolerance
            converged = true;
            return;
        end
    end
end

function Obj = CalObj(x,Zmin,RV)
    x{3}(:,x{2}) = x{1}.lower_dec(x{2});
    x= {x{1}.upper_dec,x{3}};
    Global = GLOBAL.GetObj();
    
    Zmin = min(Zmin,0);
    
    Obj = Global.problem.CalObj(x,'lower');
    Obj = max(abs(Obj-Zmin)./max(RV,eps),[],2);
end

function [c,ceq] = CalCon(x)
    x{3}(:,x{2}) = x{1}.lower_dec(x{2});
    x= {x{1}.upper_dec,x{3}};
    Global = GLOBAL.GetObj();
    c = Global.problem.CalCon(x,'lower');
    ceq = [];
end
