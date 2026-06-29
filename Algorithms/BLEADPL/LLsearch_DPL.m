function [SP,RP] = LLsearch_DPL(xu,Global)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    Nl = Global.N(2);

    l_M = Global.M(2);
    u_M = Global.M(1);
    
    
    Record = cell(1,10);
    Record_max = zeros(10,l_M);  
    Hv = zeros(1,10);
    
    if isa(xu(1),'INDIVIDUAL')
        if length(xu)<Nl
        Xl = Global.problem.Init(Nl-length(xu),'lower');
        RP =[xu,INDIVIDUAL(repmat(xu(1).udec,size(Xl,1),1),Xl,'lower')];
        else
            RP = xu;
        end
    else       
        RP = Global.Initialization(xu,Nl);
    end
    
    RP = Global.Evaluate(RP,'lower');
    
    notermination = true;
    gen = 1;
    
    [Z,~] = UniformPoint(Nl,l_M);
    
    [W,~] = UniformPoint(Nl,u_M);
    
    [FrontNo,MaxFNo_sp] = NDSort(RP.lower_objs,RP.lower_cons,Nl);
    [~,RP] = RP.adds(FrontNo');

    SP = RP;
    
    Record{1} = RP(RP.adds ==1);
    Record_max(1,:) = max(RP(RP.adds ==1).lower_objs);
    ReferencePoint0 = max(Record_max);
    
    while notermination
        %% Generate the offspring
         
        Pl = SP.lower_decs;
        I = zeros(Nl,2);
        for i= 1:Nl
           I(i,:) =  randperm(Nl,2);
        end
        P1 = Pl(I(:,1),:);
        P2 = Pl(I(:,2),:);
        Ql = DE_rand_1(Pl,P1,P2,'lower');


        SQ = INDIVIDUAL({repmat(xu,size(Ql,1),1),Ql});        
        SQ = Global.Evaluate(SQ,'lower');
        
        SR = [SP,SQ,RP(~ismember(RP.lower_decs,SP.lower_decs,'rows'))];
        
        [FrontNo,MaxFNo_rp] = NDSort(SR.lower_objs,SR.lower_cons,Nl);
        [~,SR] = adds(SR,FrontNo');
        % Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(SR.lower_objs,FrontNo);
        
        
        %% Selection
        Next_RP = FrontNo < MaxFNo_rp;                         
        % Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo_rp);


        Zmin = min(SR(all(SR.lower_cons<=0,2)).lower_objs,[],1);
        if isempty(Zmin)
            Zmin = min(SR.lower_objs,[],1);
        end
        Choose = LastSelection(SR(Next_RP).lower_objs,SR(Last).lower_objs,Nl-sum(Next_RP),Z,Zmin);
        Next_RP(Last(Choose)) = true;
        
        

       MaxFNo_sp = 1;
       while sum(FrontNo<=MaxFNo_sp)<Nl
           MaxFNo_sp = MaxFNo_sp+1;
       end
        
        
       if MaxFNo_sp>1 ||sum(FrontNo <= 1)== Nl           
           Next_SP = FrontNo < MaxFNo_sp;
           Last     = find(FrontNo==MaxFNo_sp);
           [~,Rank] = sort(CrowdDis(Last),'descend');
           Next_SP(Last(Rank(1:Nl-sum(Next_SP)))) = true;
                     
            
       else
            index =find(FrontNo==1);
            for i = 1:length(index)
               if isempty(SR(index(i)).upper_objs)
                   SR(index(i)) = Global.Evaluate(SR(index(i)),'upper');
               end
            end
            Candidate = SR(index);
            
            [u_FrontNo,u_MaxFNo] = NDSort(Candidate.upper_objs,[Candidate.upper_cons,Candidate.lower_cons],Nl);
            next = u_FrontNo<u_MaxFNo;
            
            Zmin = min(Candidate(all(Candidate.upper_cons<=0,2)).upper_objs,[],1);
            if isempty(Zmin)
                Zmin = min(Candidate.upper_objs,[],1);
            end
            last     = find(u_FrontNo==u_MaxFNo);
            Choose = LastSelection(Candidate(next).upper_objs,Candidate(last).upper_objs,Nl-sum(next),W,Zmin);
            next(last(Choose)) = true;
                    
            Next_SP = false(1,length(SR));
            Next_SP(index(next)) = true;           
       end
        RP = SR(Next_RP);
        SP = SR(Next_SP);
      
        %% termination check
        k = mod(gen,10);
        if k == 0
            k = 10;
        end
        Record{k} = RP(RP.adds==1);
        Record_max(k,:) = max(RP(RP.adds==1).lower_objs);     
        
            ReferencePoint = max(Record_max);            
            
            if sum(Hv)==0 || sum(ReferencePoint~=ReferencePoint0)>0
               ReferencePoint0 = ReferencePoint;               
               for i = 1 : min(gen,10)
                  Hv(i) = HV(Record{i}.lower_objs,ReferencePoint); 
               end
            else
                Hv(k) = HV(Record{k}.lower_objs,ReferencePoint);
            end
            
            HV_max = max(Hv(1 : min(gen,10)));
            HV_min = min(Hv(1 : min(gen,10)));      
        if gen >= 10              
            improverate = (HV_max - HV_min)/(HV_max + HV_min);
            
            notermination = improverate > 0.001;                     
        end
            
        gen = gen + 1;
            
            
    end
    
    index = find(SP.adds == 1);
    for i = 1:length(index)
        if isempty(SP(index(i)).upper_objs)
            SP(index(i)) = Global.Evaluate(SP(index(i)),'upper');
        end
    end
    
    [site,Loc] = ismember(RP.lower_decs,SP.lower_decs,'rows');
    RP(site) = SP(Loc(site));
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end