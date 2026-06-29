function [SP,RP] = LLsearch_v20(xu,Global)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    Nl = Global.N(2);
    m = Global.M(2);
    Record = cell(1,10);
    Record_max = zeros(10,m);  
    Hv = zeros(1,10);
    
    if isa(xu(1),'INDIVIDUAL')
        RP = xu(1:Nl);
        Xu = xu.upper_decs;
    else
        Xu = repmat(xu,Nl,1);   
        RP = Global.Evaluate(Global.Initialization(Xu,Nl),'lower'); %ReferencePoints
    end
    
   
    
    notermination = true;
    gen = 1;
    
    [Z,Nl] = UniformPoint(Nl,m);
    
    
    
    [FrontNo,~] = NDSort(RP.lower_objs,RP.lower_cons,inf);
    RP = RP.adds(FrontNo','lower');
    CrowdDis = CrowdingDistance(RP.lower_objs,FrontNo);
%     Next_SP = FrontNo < MaxFNo_sp;
%     Last     = find(FrontNo==MaxFNo_sp);
%     [~,Rank] = sort(CrowdDis(Last),'descend');
%     Next_SP(Last(Rank(1:Nl-sum(Next_SP)))) = true;
%     SP = RP(Next_SP);
%     FrontNo = FrontNo(Next_SP);
%     CrowdDis = CrowdDis(Next_SP);

    SP = RP;
    
    Record{1} = RP;
    Record_max(1,:) = max(RP.lower_objs);
    ReferencePoint0 = max(Record_max);
    
%     l_PF = Global.problem.lower_PF(Xu(1,:));
    while notermination
        %% Display 
%         if gen>400
%         if gen>1
%             delete(subplot(2,2,2));
%             delete(subplot(2,2,4));
%         end
%         subplot(2,2,2);        
%         title('Lower-level Obj')
%         Draw(RP.lower_objs);
% %         Draw(l_PF,'r');
%         I = sum(SP.lower_cons>0,2)<=0;
%         Draw(SP(I).lower_objs,'go');
%         Draw(SP(~I).lower_objs,'ro');        
%         
%         subplot(2,2,4);        
%         title('Lower-level Dec')
%         Draw(Xu(1,:),'b');
%         Draw(SP(I).lower_decs,'g');
%         Draw(SP(~I).lower_decs,'r');
%         
%         pause(0.0000001)
%         end
        %% Generate the offspring
         MatingPool = TournamentSelection(2,Nl,FrontNo,-CrowdDis);
         Pl = SP(MatingPool).lower_decs;
%         Pl = SP.lower_decs;
        I = zeros(Nl,2);
        for i= 1:Nl
           I(i,:) =  randperm(Nl,2);
        end
        P1 = Pl(I(:,1),:);
        P2 = Pl(I(:,2),:);
        Ql = DE_current_1_bin(Pl,P1,P2,'lower');
        
%         Ql = GA(Pl,'lower',{0.9,15,0.1,1});
%         Ql = GA(Pl,'lower');
        SQ = Global.Evaluate(INDIVIDUAL({Xu,Ql}),'lower');       
       
        SR = [SP,SQ,RP(~ismember(RP.lower_decs,SP.lower_decs,'rows'))];
        
        [FrontNo,~] = NDSort(SR.lower_objs,SR.lower_cons,inf);
        SR = SR.adds(FrontNo','lower');
        
        % Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(SR.lower_objs,FrontNo);
        
        
        %% Selection
        MaxFNo_rp = 1;
        while sum(FrontNo<=MaxFNo_rp)<Nl
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
        
        
       if MaxFNo_sp>1 ||sum(FrontNo <= 1)== Nl           
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
            
            Zmin = min(Candidate(all(Candidate.upper_cons<=0,2)).upper_objs,[],1);
            if isempty(Zmin)
                Zmin = min(Candidate.upper_objs,[],1);
            end
            last     = find(u_FrontNo==u_MaxFNo);
            Choose = LastSelection(Candidate(next).upper_objs,Candidate(last).upper_objs,Nl-sum(next),Z,Zmin);
            next(last(Choose)) = true;
            
%             u_CrowdDis = CrowdingDistance(Candidate.upper_objs,u_FrontNo);
%             last     = find(u_FrontNo==u_MaxFNo);
%             [~,rank] = sort(u_CrowdDis(last),'descend');
%             next(last(rank(1:Nl-sum(next)))) = true;
%             
            Next_SP = false(1,length(SR));
            Next_SP(index(next)) = true;           
       end
        RP = SR(Next_RP);
        SP = SR(Next_SP);
        FrontNo = FrontNo(Next_SP);
        CrowdDis = CrowdDis(Next_SP);
      
        %% termination check
        k = mod(gen,10);
        if k == 0
            k = 10;
        end
        Record{k} = RP;
        Record_max(k,:) = max(RP.lower_objs);     
        
            ReferencePoint = max(Record_max);            
            
            if gen==2 || sum(ReferencePoint~=ReferencePoint0)>0
               ReferencePoint0 = ReferencePoint;               
               for i = 1 : min(gen,10)
                  maxFNo = max(Record{i}.lower_adds);
                  Hv(i) = HV(Record{i}(Record{i}.lower_adds==maxFNo).lower_objs,ReferencePoint);
               end
            else
                maxFNo = max(Record{k}.lower_adds);
                Hv(k) = HV(Record{k}(Record{k}.lower_adds==maxFNo).lower_objs,ReferencePoint);
            end
            
            HV_max = max(Hv(1 : min(gen,10)));
            HV_min = min(Hv(1 : min(gen,10)));      
        if gen >= 10              
            improverate = (HV_max - HV_min)/(HV_max + HV_min);
            notermination = improverate > 0.001;% && gen<40*Global.l_D;                     
        end
            
        gen = gen + 1;
            
            
    end
    
    index = find(SP.lower_adds == 1);
    SP(index) = Global.Evaluate(SP(index),'upper');
    
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