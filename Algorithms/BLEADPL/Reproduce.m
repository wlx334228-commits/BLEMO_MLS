function [SQ,RQs] = Reproduce(u_Pop,Idx_ul,l_Pop,Arch_partion,Global)

%% setting
    Nu = Global.N(1);

    [S,K] = size(Arch_partion);
    
%% loading data in l_Pop
    l_xu = zeros(Nu,Global.D(1));
    for i = 1:Nu
        l_xu(i,:) = l_Pop{i}(1).upper_dec;
    end
    
%% Select the elite UL population from SQ for LL search    
    Pu = u_Pop.upper_decs;
    Pa = l_xu;
    Idx = zeros(Nu,2);
    
    K_index = 1:K;
    site =sum(ismember(Arch_partion,Idx_ul))>0; %Pop的上层向量所在分区
    if  max(K_index(site))-min(K_index(site))<=3 %Pop的上层向量所在分区数量超过2
        
        I1 = repmat(Arch_partion(1,:),S,1);
        I1 = I1(:);
        P1 = Pa(I1(:),:);
        RR = 0.5;
        
    else
        P1 = Pu;
        I1 = Idx_ul;
        RR = 0.5;  %以0.5的概率在同一分区内选择交配池
    end
    if S>1
        site = rand(Nu,2)<RR;
        for i=1:Nu
            index = 1:length(l_Pop);
            index(I1(i))=[];
            temp = index(randi(length(l_Pop)-1,1,2));
            
            k = ceil(i/S);
            index = (Arch_partion(:,k))';
            index = index(index~=I1(i));
            index = index(randi(length(index),1,2));
            temp(site(i,:))=index(site(i,:));
            Idx(i,:) = temp;
        end
    else
        for i = 1:Nu
            index = 1:length(l_Pop);
            site = ismember(Arch_partion(:),Idx_sp(i));
            index = index(site);
            Idx(i,:) = index(randi(length(index),1,2));
        end
    end
    P2 = l_xu(Idx(:,1),:);
    P3 = l_xu(Idx(:,2),:);
     
    Qu = DE_rand_1(P1,P2,P3,'upper');
                             
    I = ismember(Qu,l_xu,'rows');
    Qu(I,:) =[];
    Qu = unique(Qu,'rows');
    
    %% Lower-level Search
    if ~isempty(Qu)
        SQ = cell(1,size(Qu,1));
        RQs= cell(1,size(Qu,1));
        for i = 1:size(Qu,1)
            
            [SQ{i},RQs{i}] = LLsearch_DPL(Qu(i,:),Global);
            
        end
    else
        SQ = [];
        RQs = [];
    end
end

