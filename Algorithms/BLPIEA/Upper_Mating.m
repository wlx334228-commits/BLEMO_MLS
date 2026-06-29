function [Target,DEPool] = Upper_Mating(MatingPop,N,W,T,num)

%     N = length(MatingPop);
    K = size(W,1);
    S = N/K;
    
    inFlabels = MatingPop.labels & ~MatingPop.lower_feasibles;
    uFlabels = MatingPop.upper_feasibles;
    uFlabels(inFlabels) = false;
    Flabels = find(uFlabels);
    
    if length(Flabels)>N
        PopObj = MatingPop(Flabels).upper_objs;
        Zmin = min(PopObj,[],1);
        Zmax = max(PopObj,[],1);
        PopObj = (PopObj-repmat(Zmin,length(Flabels),1))./(repmat(max(Zmax-Zmin,1e-6),length(Flabels),1));
        [~,ZoneNum] = max(1-pdist2(PopObj,W,'cosine'),[],2);
        B = 1-pdist2(W,W,'cosine');
        [~,B] = sort(B,2,'descend');
        
        current = [];
        for i = 1:K
            current = cat(1,current,ZoneNum'==i);
        end
        current_N = sum(current,2);
        
        Target = [];
        for i = 1:K
            t=1;
            while sum(current_N(B(i,1:t)))==0
                t = t + 1;
            end
            currentID = [];
            for j = B(i,1:t)
                currentID = cat(2,currentID,find(ZoneNum'==j));
            end
            L = length(currentID);
            if L < S
                currentID = [currentID,currentID(randi(L,1,S-L))];
            else
                [~,IndexP] = NSGAII_Update(MatingPop(Flabels(currentID)),S,'upper');
                currentID = currentID(IndexP);
            end
            
            Target = cat(2,Target,Flabels(currentID));
            
            ZoneNum(currentID) = 0;
            current(:,currentID) = 0;
            current_N = sum(current,2);
            
        end
    else
        Target = Flabels;
        FrontNo = NDSort(MatingPop.upper_objs,MatingPop.upper_cons,inf);
        CrowdDis = CrowdingDistance(MatingPop.upper_objs,FrontNo);
        IndexP = TournamentSelection(2,N-length(Flabels),FrontNo,-CrowdDis);
        Target = cat(2,Target,IndexP);
    end
    
    B = pdist2(MatingPop(Target).upper_objs,MatingPop.upper_objs);
    B(B==0) = inf;
    [~,B] = sort(B,2);
    B = B(:,1:T);
    DEPool = [];
    for i=1:N
        if rand < 0.5
            DEPool = cat(1,DEPool,B(i,randperm(T,num)));
        else
            DEPool = cat(1,DEPool,randperm(length(MatingPop),num));
        end
    end
    
    Target = Target';
end