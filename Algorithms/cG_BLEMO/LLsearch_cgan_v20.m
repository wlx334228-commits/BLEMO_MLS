function [SP,RP] = LLsearch_cgan_v20(N,rp1,cgan,Global)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    
    Nl = Global.N(2);
    xu = rp1(1).upper_dec;
    m  = Global.M(2);       
    [Z,Nl] = UniformPoint(Nl,m);
    l_PF = Global.problem.lower_PF(xu);
    
    
%     No =10;
%     N  =No;
    %% 初始化

     Ux = repmat(xu,Nl-length(rp1),1);
     RP = [rp1,Global.Evaluate(INDIVIDUAL({Ux,cgan.Generator(Ux)}),'lower')];

    
    [FrontNo,~] = NDSort(RP.lower_objs,RP.lower_cons,inf);  
    RP = RP.adds(FrontNo','lower');
    
%     Next_SP = FrontNo < MaxFNo_sp;
%     % Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(RP.lower_objs,FrontNo);
%     % Select the solutions in the last front based on their crowding distances
%     Last     = find(FrontNo==MaxFNo_sp);
%     [~,Rank] = sort(CrowdDis(Last),'descend');
%     Next_SP(Last(Rank(1:Nl-sum(Next_SP)))) = true;  
%     SP = RP(Next_SP);
%     
%     FrontNo = FrontNo(Next_SP);
%     CrowdDis = CrowdDis(Next_SP);
    Next_SP= true(1,Nl);     
    SP = RP;

   %终止条件初始化 
    notermination = true;    
    gen = 1;
    
    Record = cell(1,10);
    Record_max = zeros(10,m);
    Record{1} = RP;
    Record_max(1,:) = max(RP.lower_objs);
    Hv = zeros(1,10);
    FNos = zeros(Nl,10);
    FNos(:,1) = RP.lower_adds;
    
%     l_PF = Global.problem.lower_PF(SP(1).upper_dec);
    while notermination
               %% Display 
%         if gen>1
%             delete(subplot(2,2,2));
%             delete(subplot(2,2,4));
%         end
%         subplot(2,2,2);        
%         title('Lower-level Obj')
%         Draw(l_PF,'r');
%         Draw(RP.lower_objs);
%         I = SP.lower_cons<=0;
%         Draw(SP(I).lower_objs,'go');
%         Draw(SP(~I).lower_objs,'ro');        
%         
%         subplot(2,2,4);        
%         title('Lower-level Dec')
%         Draw(SP(1).upper_dec,'b');
%         Draw(SP(I).lower_decs,'g');
%         Draw(SP(~I).lower_decs,'r');
%         
%         pause(0.0000001)


        % Update matepool      
         MatingPool = TournamentSelection(2,N,FrontNo,-CrowdDis);
         Pl = SP(MatingPool).lower_decs;
        %% Generate the offspring         
            idx = zeros(N,2);
            for i = 1:N
               idx(i,:) = randperm(Nl,2);
            end
            P1 = SP(idx(:,1)).lower_decs;
            P2 = SP(idx(:,2)).lower_decs;
            Q = DE_current_1_bin(Pl,P1,P2,'lower');


        SQ = Global.Evaluate(INDIVIDUAL({repmat(xu,size(Q,1),1),Q}),'lower');

        SR = [SQ,SP,RP(~ismember(RP.lower_decs,SP.lower_decs,'rows'))];
        
        %% Non-dominated sorting
        [FrontNo,~] = NDSort(SR.lower_objs,SR.lower_cons,inf);
        SR = SR.adds(FrontNo','lower');
        
        %% Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(SR.lower_objs,FrontNo);
        
        %% Selection
        MaxFNo_rp = 1;
        while sum(FrontNo <= MaxFNo_rp)<Nl
            MaxFNo_rp = MaxFNo_rp+1;
        end
        Next_RP = FrontNo < MaxFNo_rp;                         
        % Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo_rp);
%         [~,Rank] = sort(CrowdDis(Last),'descend');
%         Next_RP(Last(Rank(1:Nl-sum(Next_RP)))) = true;

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
        if MaxFNo_sp>1 ||sum(FrontNo <= MaxFNo_sp)== Nl           
           Next_SP = FrontNo < MaxFNo_sp;
           Last     = find(FrontNo==MaxFNo_sp);
           [~,Rank] = sort(CrowdDis(Last),'descend');
           Next_SP(Last(Rank(1:Nl-sum(Next_SP)))) = true;
           else
            index =find(FrontNo==1);
            SR(index) = Global.Evaluate(SR(index),'upper');
            Candidate = SR(index);
            
            [u_FrontNo,u_MaxFNo] = NDSort(Candidate.upper_objs,Candidate.upper_cons,Nl);
            next = u_FrontNo<u_MaxFNo;
%             u_CrowdDis = CrowdingDistance(Candidate.upper_objs,u_FrontNo);
%             last     = find(u_FrontNo==u_MaxFNo);
%             [~,rank] = sort(u_CrowdDis(last),'descend');
%             next(last(rank(1:Nl-sum(next)))) = true;

            Zmin = min(Candidate(all(Candidate.upper_cons<=0,2)).upper_objs,[],1);
            if isempty(Zmin)
                Zmin = min(Candidate.upper_objs,[],1);
            end
            last     = find(u_FrontNo==u_MaxFNo);
            Choose = LastSelection(Candidate(next).upper_objs,Candidate(last).upper_objs,Nl-sum(next),Z,Zmin);
            next(last(Choose)) = true;
            
            Next_SP = false(1,length(SR));
            Next_SP(index(next)) = true;           
        end
        %% Population for next generation
        RP = SR(Next_RP);
        SP = SR(Next_SP);
        FrontNo = FrontNo(Next_SP);
        CrowdDis = CrowdDis(Next_SP);
        
        gen = gen + 1;   
         %% termination check
         % Record the last 10 generations of populations
        k = mod(gen,10);
        if k == 0
            k = 10;
        end
        Record{k} = RP;
        Record_max(k,:) = max(RP.lower_objs); 
        FNos(:,k) = RP.lower_adds;
         
        if gen >= 10+max(0,10-N)
            ReferencePoint = max(Record_max); 
            if sum(Hv)==0|| sum(ReferencePoint~=ReferencePoint0)>0
               ReferencePoint0 = ReferencePoint;
               for i = 1 : min(gen,10)
                  maxFNo = max(FNos(:,i));
                  Hv(i) = HV(Record{i}(FNos(:,i)==maxFNo).lower_objs,ReferencePoint);
               end
            else
                maxFNo = max(FNos(:,k));
                Hv(k) = HV(Record{k}(FNos(:,k)==maxFNo).lower_objs,ReferencePoint);
            end
            
            HV_max = max(Hv);
            HV_min = min(Hv);
            
            
            notermination = (HV_max - HV_min)/(HV_max + HV_min) > 0.001 ;
%             N =max(N,No + max(0,ceil(-log10(max((HV_max - HV_min)/(HV_max + HV_min),0.001)))*5));
        end
                      
    end
    index = find(FrontNo == 1);
    SP(index) = Global.Evaluate(SP(index),'upper');
    
    [site,Loc] = ismember(RP.lower_decs,SP.lower_decs,'row');
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