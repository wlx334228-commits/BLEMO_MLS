function [output,Idx_oA,Archive,RPs] = RefineSearch(output,Idx_oA,Archive,RPs,Global)


        Nu = Global.N(1);
        Nl = Global.N(2);

        
        [Z,~] = UniformPoint(Nl,Global.M(1));
       
        %% 

    for t = 1:20    
        isUpdate = false;


        Au = zeros(length(Archive),Global.D(1));        
        for i = 1:length(Archive)  
            Au(i,:) = Archive{i}(1).upper_dec;
            
            idx_p = find(Idx_oA == i);           
            n = length(idx_p);       
            xu = output(idx_p(1)).upper_dec;

            p1 = output(idx_p).lower_decs;
            site = zeros(n,2);
            for j = 1:n
               site(j,:) = randperm(Nl,2);
            end
            p2 = Archive{i}(site(:,1)).lower_decs;           
            p3 = Archive{i}(site(:,2)).lower_decs;
            
            Xl = DE_rand_1(p1,p2,p3,'lower'); 
  
            off = INDIVIDUAL({repmat(xu,size(Xl,1),1),Xl}); 
            off = Global.Evaluate(off,'lower');
            %
            r = [Archive{i},RPs{i}(~ismember(RPs{i}.lower_decs,Archive{i}.lower_decs,'rows')),off];
            
            % Non-dominated sorting
            [l_FrontNo,maxFNo_rp] = NDSort(r.lower_objs,r.lower_cons,Nl);
            [~,r] = adds(r,l_FrontNo');   
            
            off_FN = l_FrontNo(length(r)-length(off)+1:end);             
                       
            if sum(off_FN == 1)>0               
                off(off_FN == 1) = Global.Evaluate(off(off_FN == 1),'upper');
                output = [output,off(off_FN == 1)];
                Idx_oA = [Idx_oA;i*ones(sum(off_FN == 1),1)];
                isUpdate = true;
            end 
            r(length(r)-length(off)+1:end) = off;
            
            Next_RP = l_FrontNo < maxFNo_rp;
            l_CrowdDis = CrowdingDistance(r.lower_objs,l_FrontNo);
            Last     = find(l_FrontNo==maxFNo_rp);
            [~,Rank] = sort(l_CrowdDis(Last),'descend');
            Next_RP(Last(Rank(1:Nl-sum(Next_RP)))) = true;
            RPs{i} = r(Next_RP);
            
            maxFNo_sp = 1;
            while sum(l_FrontNo<=maxFNo_sp)<Nl
                maxFNo_sp = maxFNo_sp+1;            
            end
            if maxFNo_sp == 1
                Archive{i} = r(l_FrontNo == 1);              
            else
               Archive{i} = r(Next_RP);
            end
            clear r off             
        end
        
        %% Update output & Archive
        if isUpdate
            [FrontNo,~] = NDSort(output.upper_objs,[output.upper_cons,output.lower_cons],1);
            output = output(FrontNo==1);
            if length(output)>Global.N(3)
                Zmin = min(output.upper_objs);
                [W,~]= UniformPoint(Global.N(3),2);
                [~,I] = max(1-pdist2(output.upper_objs-Zmin,W,'cosine'));
                I = unique(I);
                output = output(I);
            end
            A = Archive;
            U = RPs;
            o_udecs = output.upper_decs;
    
            [~,Loc] = ismember(o_udecs,Au,'rows');
            I = unique(Loc);
            Archive = cell(1,length(I));
            RPs = cell(1,length(I));
            for i = 1:length(I)
                Archive{i} = A{I(i)};
                RPs{i} = U{I(i)};
            end
            [~,Idx_oA] = ismember(Loc,I);
            
         for i = 1:length(Archive)
             if length(Archive{i}) > Nl
                 [fn,MaxFNo] = NDSort(Archive{i}.upper_objs,Archive{i}.upper_cons,Nl);
                 zmin = min(Archive{i}(all(Archive{i}.upper_cons<=0,2)).upper_objs,[],1);
                 if isempty(zmin)
                     zmin = min(Archive{i}.upper_objs,[],1);
                 end
                 
                 Next = fn<MaxFNo;
                 Last     = find(fn==MaxFNo);
                 Choose = LastSelection(Archive{i}(Next).upper_objs,Archive{i}(Last).upper_objs,Nl-sum(Next),Z,zmin);
                 Next(Last(Choose)) = true;
                 
                 Archive{i} = Archive{i}(Next);
             end
         end
        end
    end
end


function Choose = LastSelection(outputObj1,outputObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    outputObj = [outputObj1;outputObj2] - repmat(Zmin,size(outputObj1,1)+size(outputObj2,1),1);
    [N,M]  = size(outputObj);
    N1     = size(outputObj1,1);
    N2     = size(outputObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(outputObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = outputObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(outputObj,[],1)';
    end
    % Normalization
    outputObj = outputObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(outputObj,Z,'cosine');
    Distance = repmat(sqrt(sum(outputObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
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