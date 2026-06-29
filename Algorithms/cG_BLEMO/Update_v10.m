function [Pop,FNo_p,Idx_sp,Archive,Arch_partion] = Update_v10(R,Idx_sr,SR,W,S)

    [K,u_M] = size(W);
  
    
   %% Update Pop
    % Allocation of solutions to subproblems
    % Transformation
%     P = R.upper_objs-min(R.upper_objs);
%     P = P./repmat(sqrt(sum(P.^2,2)),1,2);
%     
%     Zmin = min(P);
%     [ad,transformation] = max(1-pdist2(P-Zmin,W,'cosine'),[],2);

[ad,transformation] = max(1-pdist2(R.upper_objs-min(R.upper_objs),W,'cosine'),[],2);
    
    partition = zeros(S,K);
    for i = 1 : K
        if i<K
            [z,S] = UniformPoint(S+1,u_M);  
            S = S - 1;
            z = z(1:S,:);  
        else
            [z,S] = UniformPoint(S,u_M);  
        end
        
        current = find(transformation==i);
        site = ~ismember(current,partition);
        current  = current(site);
        
        if length(current) < S
            current = [current;zeros(S-length(current),1)];
        else
            r = R(current);
            [FNo_current,maxFNo] = NDSort(r.upper_objs,[r.upper_cons,r.lower_cons],S);
           
            Last = find(FNo_current==maxFNo);
            
            % NSGA III
            zmin = min(r(all(r.upper_cons<=0,2)).upper_objs,[],1);
            if isempty(zmin)
                zmin = min(r.upper_objs,[],1);
            end
            
            Next = FNo_current < maxFNo;
            Choose = LastSelection(r(Next).upper_objs,r(Last).upper_objs,S-sum(Next),z,zmin);
            Next(Last(Choose)) = true;
            
            current = current(Next);
            
        end
        partition(:,i) = current;
    end
    
    index = (1:length(R))';
    site = ~ismember(index,partition);
    candidate = index(site);
    site = partition==0;
    partition(site) = candidate(randperm(length(candidate),min(length(candidate),sum(site(:)))));     
    
    Pop = R(partition(:));
    [FNo_p,~] = NDSort(Pop.upper_objs,[Pop.upper_cons,Pop.lower_cons],inf);
    Idx_temp = Idx_sr(partition);
    if size(partition,1)==1
        Idx_temp = Idx_temp';
    end
    
    Idx_tp = unique(Idx_temp(:));      
    site = ismember(Idx_sr,Idx_tp);
    R(site) = [];
    Idx_sr(site) = [];
    transformation(site) = [];
    ad(site) = [];

    
    %% Update Archive
    partition = zeros(S,K);
    
    for i = 1 : K             
        Idx_selected = Idx_temp(:,i);
        site = ~ismember(Idx_selected,partition);
        Idx_selected = unique(Idx_selected(site),'stable');
        
        current = transformation==i;
        Candidate = R(current);
        Idx_candidate = Idx_sr(current);
        ad_candidate  = ad(current);     
        
        Selected = Pop((i-1)*S+1:i*S);
        
        
        if ~isempty(Idx_candidate)
            r = [Selected,Candidate];
            zmin = min(r.upper_objs);
            [FNo_candidate,~] = NDSort(Candidate.upper_objs,[Candidate.upper_cons,Candidate.lower_cons],inf);
        end
        while length(Idx_selected) < S          
            if isempty(Idx_candidate)
                % Randomly select solutions and join to the current subproblem
                Idx = 0;
            else                           
                [minFNo_c,~] = min(FNo_candidate);
                site = find(FNo_candidate == minFNo_c);
                if length(site)>1
                    [ad_S,~] = max(1-pdist2(Candidate(site).upper_objs-zmin,Selected.upper_objs-zmin,'cosine'),[],2);
                    [~,I1] = min(ad_S);
                    I = site(I1);
                else
                    I = site;
                end
                Idx = Idx_candidate(I);
                
                site = Idx_candidate == Idx;
                Selected = [Selected,Candidate(site)];
                Candidate(site) = [];
                Idx_candidate(site) = [];
                FNo_candidate(site) = [];
                ad_candidate(site) = [];

            end
            site = Idx_sr == Idx;
            R(site) = [];
            Idx_sr(site) = [];
            transformation(site) = [];
            ad(site) = [];
            Idx_selected =  [Idx_selected;Idx];
        end
        partition(:,i) = Idx_selected;
    end
    
    Idx_add = unique(Idx_sr); 
    site = partition ==0;
    partition(site) = Idx_add(randperm(length(Idx_add),sum(site(:))));
    
    Archive_next = partition(:);
     
    [~,Arch_partion] = ismember(partition,Archive_next);    
    [~,Idx_sp] = ismember(Idx_temp(:),Archive_next);
    
    na = K*S;
    Archive = cell(1,na);

    for i = 1:na
        Archive{i} = SR{Archive_next(i)};

    end
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