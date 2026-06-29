function [SQ,RQs] = Reproduce_lr(N,tmax,u_Pop,Idx_ul,l_Pop,Arch_partion,cgan,W,Global)
%UNTITLED3 姝ゅ鏄剧ず鏈夊叧姝ゅ嚱鏁扮殑鎽樿
%   姝ゅ鏄剧ず璇︾粏璇存槑

%% setting
    Nu = Global.Nu;
    Nl = Global.Nl;

    [S,K] = size(Arch_partion);
    
    [Z,S] = UniformPoint(S,Global.u_M); 
    
    r =[];
%% loading data in l_Pop
    l_xu = zeros(Nu,Global.u_D);
    for i = 1:Nu
        l_xu(i,:) = l_Pop{i}(1).udec;
    end
    
%% Select the elite UL population from SQ for LL search    
    pop = u_Pop;
    Pa = l_xu;
    
    for t=1:tmax
        Pu = pop.udecs;  
        
        K_index = 1:K;
        site =sum(ismember(Arch_partion,Idx_ul))>0; %Pop的上层向量所在分区
        if  max(K_index(site))-min(K_index(site))<=3 %Pop的上层向量所在分区数量超过2
            
            I1 = repmat(Arch_partion(1,:),S,1);
            I1 = I1(:);
            P1 = Pa(I1(:),:);
            RR = 0.7;
            
        else
            P1 = Pu;
            
            RR = 0.7;  %以0.7的概率在同一分区内选择交配池
        end
        
        Idx = zeros(K*S,2);
        for i = 1:K*S
            if rand<RR
                Index = Arch_partion(:,ceil(Idx_ul(i)/S));
            else
                Index = Arch_partion(:);
            end
            Index = Index(Index~=Idx_ul(i));
            Idx(i,:) = Index(randi(length(Index),1,2));
        end
        P2 = Pa(Idx(:,1),:);
        P3 = Pa(Idx(:,2),:);
        
        Qu = DE(P1,P2,P3,'upper');
        
        
        
        Xl = cgan.Generator(Qu);
        off = INDIVIDUAL(Qu,Xl,'upper');
        R = [pop,off];
        [pop,Idx_ul,Pa,Arch_partion] = m2mSelection(R,W,Z,Nu,K,S);
        
        r = [r,off];
        if length(r)>Nu
            [r,~,~,~] = m2mSelection(r,W,Z,Nu,K,S);
        end              
        
        %% Display
%         if t>1
%             delete(subplot(2,2,1));
%             delete(subplot(2,2,3));
%         end
%         
%         subplot(2,2,1);
%         title('Upper-level Obj')
%         
%         Draw(u_Pop.uobjs,'go');
%         Draw(pop.uobjs,'yo');
%         
%         subplot(2,2,3);
%         title('Upper-level Dec')
%         Draw(l_xu);
%         Draw(pop.udecs,'y');
%         Draw(u_Pop.udecs,'g');
%         
%         hold off;
%         pause(0.00001);      
        
    end
        
        zmin = min(r.uobjs);
        [~,transformation] = max(1-pdist2(r.uobjs-zmin,W,'cosine'),[],2);
    
        R = [u_Pop,r];
        [FrontNo,~] = NDSort(R.uobjs,R.ucons,inf);
        r_FNo = FrontNo(length(u_Pop)+1:end);
        r_CrowdDis = CrowdingDistance(r.uobjs,r_FNo);             
        
        site = r_FNo == 1;
        elite = r(site);
        st = transformation(site);
        
        ct = transformation(~site);
        r = r(~site);
        r_FNo = r_FNo(~site);
        r_CrowdDis = r_CrowdDis(~site);
        if length(elite)<ceil(Nu/2)
            for i = 1:K
                site =  ct == i;
                if sum(site)>0&&sum(st==i)<Nu/(2*K)
                    current = r(site);
                    c_FNo = r_FNo(site);
                    c_CrowdDis = r_CrowdDis(site);
                    while sum(st==i)<Nu/(2*K)
                        if length(current)>Nu/(2*K)-sum(st==i)
                            I = TournamentSelection(2,1,c_FNo,-c_CrowdDis);
                            elite = [elite,current(I)];
                            st = [st;i];
                        
                            current(I) = [];
                            c_FNo(I) = [];
                            c_CrowdDis(I)=[];
                        else
                            elite = [elite,current];
                            st = [st;ones(length(current),1)*i];
                        end
                    end
                else
                    continue;
                end
            end
           
        end

%         elite = r;
        
        Qu = elite.udecs;
    
        I = ismember(Qu,l_xu,'rows');
        Qu(I,:) =[];
        Qu = unique(Qu,'rows');
        
        %% 准确搜索下层最优解
        if ~isempty(Qu)
            SQ = cell(1,size(Qu,1));
            RQs= cell(1,size(Qu,1));
            for i = 1:size(Qu,1)
%                 if tmax>1
%                     site = ismember(elite.udecs,Qu(i,:),'rows');
%                     SQ{i} = INDIVIDUAL(elite(site),'lower'); 
%                 else
%                     SQ{i} =  INDIVIDUAL(Qu(i,:),cgan.Generator(Qu(i,:)),'lower');
%                 end
                [SQ{i},RQs{i}] = LLsearch_v20(Qu(i,:),Global);

            end
        else
            SQ = [];
            RQs = [];
        end
end

function [pop,Idx_sp,Pa,Arch_partion]= m2mSelection(R,W,Z,Nu,K,S)

    Ru = unique(R.udecs,'rows');
    [~,Idx_R] = ismember(R.udecs,Ru,'rows');
    %% Assosiation
    
    [FrontNo,MaxFNo] = NDSort(R.uobjs,R.ucons,Nu);
    P = R.uobjs-min(R.uobjs);
    P = P./repmat(sqrt(sum(P.^2,2)),1,2);
    
    Zmin = min(P);
    [~,transformation] = max(1-pdist2(P-Zmin,W,'cosine'),[],2);
    partition = zeros(S,K);
    index = (1:length(R))';
    
    %% Selection
    for i = 1 : K
        current = find(transformation==i);
        site = ~ismember(current,partition);
        current = current(site);

        if length(current) < S
            % Randomly select solutions and join to the current subproblem
            
            current = [current;zeros(S-length(current),1)];
        elseif length(current) > S
            r = R(current);
            FNo_current = FrontNo(current);
            FNo_c  = unique(FNo_current);
            j = 1;
            while sum(FNo_current<=FNo_c(j))<S
                j = j+1;
            end
            maxFNo = FNo_c(j);

            Next = FNo_current < maxFNo;
            Last = find(FNo_current==maxFNo);
            zmin = min(r(all(r.ucons<=0,2)).uobjs,[],1);
            if isempty(zmin)
                zmin = min(r.uobjs,[],1);
            end

            Choose = LastSelection(r(Next).uobjs,r(Last).uobjs,S-sum(Next),Z,zmin);
            Next(Last(Choose)) = true;
            current = current(Next);
        end
        partition(:,i) = current;
    end
    site = ~ismember(index,partition);
    idx_candidate = index(site);
    site = partition == 0;
    partition(site) = idx_candidate(randperm(length(idx_candidate),sum(site(:))));
    pop = R(partition(:));
    Temp = Idx_R(partition);
    Idx_sp = Temp(:);
    
    Arch_partion = zeros(S,K);
    for i = 1:K
        selected = Temp(:,i);
        site = ~ismember(selected,Arch_partion);
        current = unique(selected(site),'stable');
        
        No = find(transformation==i);
        Idx_candidate = Idx_R(No);
        site = ~ismember(Idx_candidate,[Temp,Arch_partion]);
        Idx_candidate = Idx_R(site);
        No = No(site);
        Candidate = R(No);
        FNo_candidate = FrontNo(No);
        Selected = pop((i-1)*S+1:i*S);
        
        while length(current)<S
            r = [Selected,Candidate];
            zmin = min(r.uobjs);
            if isempty(Candidate)
                current = [current;zeros(S-length(current),1)];
            else
                
                [minFNo,~] = min(FNo_candidate);
                I = find(FNo_candidate == minFNo);
                if length(I)>1
                    [ad,~]= max(1-pdist2(Candidate(I).uobjs-zmin,Selected.uobjs-zmin,'cosine'),[],2);
                    [~,site] = min(ad);
                    I = I(site);
                end

                current = [current;Idx_candidate(I)];
                Selected =[Selected,Candidate(I)];
                site = Idx_candidate == Idx_candidate(I);
                Candidate(site) = [];
                Idx_candidate(site) = [];
                FNo_candidate(site) = [];
            end
        end
        Arch_partion(:,i)=current;
    end
    
    Index = 1:size(Ru,1);
    site = ~ismember(Index,Arch_partion);
    if sum(site<1)
       site = true(1,length(Index));
    end
    Index = Index(site);
    
    site = Arch_partion==0;
    if sum(site(:)) <= length(Index)
        I = randperm(length(Index),sum(site(:)));
    else
        I =randi(length(Index),sum(site(:)),1);
    end
    Arch_partion(site) = Index(I);
    
    Pa = Ru(Arch_partion(:),:);
    [~,Idx_sp] = ismember(Idx_sp,Arch_partion(:));
    [~,Arch_partion] = ismember(Arch_partion,Arch_partion(:));
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
