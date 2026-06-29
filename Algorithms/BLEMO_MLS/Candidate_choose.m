function [ID,rnd] = Candidate_choose(Population,N,W)
    Global = GLOBAL.GetObj();
    
    ID = [];
    labels = Population.labels;
    ind1 = find(labels & Population.bilevel_feasibles);
    ind2 = find(~labels);
    %ind2 = find(~labels & Population.upper_feasibles);
%     ind1 = [];
    if ~isempty(ind2)
        if ~isempty(ind1)
            [I_epsilon,rate] = C_I_epsilon(Population(ind2).upper_objs,Population(ind1).upper_objs);
            ID = ind2(I_epsilon>=1|rate>=0.8*length(ind1));
        else
            ID = cat(2,ID,ind2);
        end
        
%         if ~isempty(ind1)
%             FronNo = NDSort(Population([ind2,ind1]).upper_objs,inf);
%             ID = cat(2,ID,ind2((FronNo(1:length(ind2))<=min(FronNo(1+length(ind2):end)))));
%             ind2_temp = ind2((FronNo(1:length(ind2))>min(FronNo(1+length(ind2):end))));
%             
%             if ~isempty(ind2_temp)
%                 I_epsilon = Cal_I_epsilon(Population(ind2_temp).upper_objs,Population(ind1).upper_objs);
%                 ID = cat(2,ID,ind2_temp(I_epsilon'>0.8*length(ind1)));
%             end
%         else
%             ID = cat(2,ID,ind2);
%         end
        
%         if ~isempty(ind1)
%             Zmax = max(Population(ind1).upper_objs,[],1);
%             Zmin = min(Population(ind1).upper_objs,[],1);
%             
%             M = size(Zmax,2);
%             
%             Upper_objs_ind2 = Population(ind2).upper_objs;
%             ind = all(Upper_objs_ind2<=Zmax,2)&all(Upper_objs_ind2>=Zmin,2);
%             Zone_1 = ind2(ind);
%             Zone_temp = ind2(~ind);
%             
%             if ~isempty(Zone_1)
%                 I_epsilon = Cal_I_epsilon(Population(Zone_1).upper_objs,Population(ind1).upper_objs);
%                 FronNo = NDSort(Population([Zone_1,ind1]).upper_objs,inf);
%                 ID = cat(2,ID,Zone_1((I_epsilon'>0.8*length(ind1))|(FronNo(1:length(Zone_1))...
%                     <=min(FronNo(1+length(Zone_1):end)))));
%             end
%             
%             if ~isempty(Zone_temp)
%                 Upper_objs_Zone_temp = Population(Zone_temp).upper_objs;
%                 inds = false(1,length(Zone_temp));
%                 for i = 1:M
%                     index = 1:M;
%                     index(i) = [];
%                     inds(all(Upper_objs_Zone_temp(:,index)< Zmin(index),2)) = true;
% %                     ID = cat(2,ID,Zone_temp(all(Upper_objs_Zone_temp(:,index)< Zmin(index),2)));
% %                     ID = cat(2,ID,Zone_temp((Upper_objs_Zone_temp(:,i)> Zmax(i)) & ...
% %                         all(Upper_objs_Zone_temp(:,index)< Zmin(index),2)));
%                 end
%                 ID = cat(2,ID,Zone_temp(inds));
%             end
%         else
%             ID = ind2;
%         end

    else
        ind2 = find(~labels);
        if ~isempty(ind2)
            [~,rank] = sort(Population(ind2).upper_CVs,'ascend');
            ID = ind2(rank(1));
        end
    end
    
    ID = unique(ID);
    rnd = length(ID)/max(sum(~labels),eps);
    
    if length(ID)>N
        [~,Next]=Upper_Select(Population(ID),N,W,'upper');
        ID = ID(Next);
    end
end

function [I_epsilon,rate] = C_I_epsilon(PopObj1,PopObj2)

    N1 = size(PopObj1,1);
    N2 = size(PopObj2,1);
    Min_obj = min([PopObj1;PopObj2],[],1);
    
    PopObj1 = PopObj1 - Min_obj;
    PopObj2 = PopObj2 - Min_obj;
    
    I_epsilon = zeros(1,N1);
    rate = zeros(1,N1);
    for i=1:N1
        I_epsilon(i) = min(max(PopObj2./max(repmat(PopObj1(i,:),N2,1),eps),[],2));
%         I_epsilon(i) = min(max(repmat(PopObj1(i,:),N2,1)./max(PopObj2,eps),[],2));
        rate(i) = sum(max(repmat(PopObj1(i,:),N2,1)./max(PopObj2,eps),[],2)<1);
    end
    
end