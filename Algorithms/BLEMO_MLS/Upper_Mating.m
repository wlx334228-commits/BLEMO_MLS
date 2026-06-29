function [Target,DEPool] = Upper_Mating(Pop,N,W,T,num)

    K = size(W,1);
    S = N/K;
    
%     [I,rate] = Candidate_choose(MatingPop,LM,W);
%     
    inFlabels = Pop.labels & ~Pop.lower_feasibles;
    uFlabels = Pop.upper_feasibles;
    uFlabels(inFlabels) = false;
    Flabels = find(uFlabels);
    
    if length(Flabels)>N
        MatingPop = Pop(Flabels);
        LM = length(MatingPop);
        PopObj = Pop(Flabels).upper_objs;
        Zmin = min(PopObj,[],1);
        Zmax = max(PopObj,[],1);
        PopObj = (PopObj-repmat(Zmin,LM,1))./(repmat(max(Zmax-Zmin,1e-6),LM,1));
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
            if L == S
                Target = cat(2,Target,currentID);
            elseif L < S
                Target = cat(2,Target,[currentID,currentID(randi(L,1,S-L))]);
            else
                
                FrontNo = NDSort(MatingPop(currentID).upper_objs,MatingPop(currentID).upper_cons,inf);
                CrowdDis = CrowdingDistance(MatingPop(currentID).upper_objs,FrontNo);
                
                %             if isempty(I) %|| rand < 1-rate
                IndexP = currentID(TournamentSelection(2,S,FrontNo,-CrowdDis));
                %             else
                %                 IndexP = [];
                %                 for j = 1:S
                %                     IndexP_j = randperm(L,2);
                %                     a = ismember(currentID(IndexP_j),I);
                %                     if any(a)
                %                         if all(a)
                %                             Dis = min(pdist2(MatingPop(currentID(IndexP_j)).upper_objs,MatingPop(labels).upper_objs),[],2);
                %                             [~,Ind_j] = min(Dis);
                %                             IndexP_j = currentID(IndexP_j(Ind_j));
                %                         else
                %                             Dis_a = min(pdist2(MatingPop(currentID(IndexP_j(a))).upper_objs,MatingPop(labels).upper_objs),[],2);
                %                             Dis_b = min(pdist2(MatingPop(currentID(IndexP_j(~a))).upper_objs,MatingPop(I).upper_objs),[],2);
                %                             if Dis_b<=Dis_a
                %                                 IndexP_j = currentID(IndexP_j(~a));
                %                             else
                %                                 IndexP_j = currentID(IndexP_j(a));
                %                             end
                %                         end
                %                     else
                %                         if FrontNo(IndexP_j(1))==FrontNo(IndexP_j(2))
                %                             [~,Ind_j] = min(CrowdDis(IndexP_j));
                %                             IndexP_j = currentID(IndexP_j(Ind_j));
                %                         else
                %                             [~,Ind_j] = min(FrontNo(IndexP_j));
                %                             IndexP_j = currentID(IndexP_j(Ind_j));
                %                         end
                %                     end
                %
                %                     IndexP = cat(2,IndexP,IndexP_j);
                %                 end
                %             end
                Target = cat(2,Target,IndexP);
            end
            %             [~,IndexP] = NSGAII_Update(MatingPop(currentID),S,'upper');
            
            ZoneNum(Target) = 0;
            current(:,Target) = 0;
            current_N = sum(current,2);
            
        end
        
        Target = Flabels(Target);
        
    else
        Target = Flabels;
        FrontNo = NDSort(Pop.upper_objs,Pop.upper_cons,inf);
        CrowdDis = CrowdingDistance(Pop.upper_objs,FrontNo);
        IndexP = TournamentSelection(2,N-length(Flabels),FrontNo,-CrowdDis);
        Target = cat(2,Target,IndexP);
    end
%     B = pdist2(MatingPop(Target).upper_objs,MatingPop.upper_objs);
%     B(B==0) = inf;
%     [~,B] = sort(B,2);
%     B = B(:,1:T);
%     
%     DEPool = [];
%     for i=1:N
%         if rand < 0.5
%             DEPool = cat(1,DEPool,B(i,randperm(T,num)));
%         else
%             DEPool = cat(1,DEPool,randperm(N,num));
%         end
%     end
    Target = transpose(Pop(Target));
    
end