function [redundant_dec_uobj, redundant_dec_lobj] = group(Global)

    [groups_u,groupNum_u] = group2(Global,'upper');
    [groups_l,groupNum_l] = group2(Global,'lower');
    
    redundant_dec_uobj = [];
    for i=1:groupNum_u
        if length(groups_u{i})==1 && length(groups_u{i}{:})==1
            redundant_dec_uobj=cat(2,redundant_dec_uobj,[groups_u{i}{:}]);
        end
    end
    
    if sum(redundant_dec_uobj<=Global.D(1))==Global.D(1)
        redundant_dec_uobj(redundant_dec_uobj<=Global.D(1))=[];
    end
    
    if sum(redundant_dec_uobj>Global.D(1))==Global.D(2)
        redundant_dec_uobj(redundant_dec_uobj>Global.D(1))=[];
    end
    
    redundant_dec_lobj = [];
    for i=1:groupNum_l
        if length(groups_l{i})==1 && length(groups_u{i}{:})==1
            redundant_dec_lobj=cat(2,redundant_dec_lobj,[groups_l{i}{:}]);
        end
    end
    
    if sum(redundant_dec_lobj<=Global.D(1))==Global.D(1)
        redundant_dec_lobj(redundant_dec_uobj<=Global.D(1))=[];
    end
    
    if sum(redundant_dec_lobj>Global.D(1))==Global.D(2)
        redundant_dec_lobj(redundant_dec_lobj>Global.D(1))=[];
    end
    
end

function [groups,groupNum] = group2(Global,type)

epsilon = 1e-4;

xrange = [Global.upper_domain,Global.lower_domain];
dim = sum(Global.D);
u_dim = Global.D(1);

xL = xrange(1,:)+rand(1,dim).*(xrange(2,:)-xrange(1,:));

x1 = xrange(1,:)+1e-1;
Pop1 = Global.Evaluate(INDIVIDUAL({x1(1:u_dim),x1(u_dim+1:dim)}),type);
fit1 =[Pop1.([type,'_obj']),Pop1.([type,'_con'])];

num = size(fit1,2);
archiveL = -1*ones(dim,num);
lambdaObjL = -1*ones(dim,dim,Global.M(2));
lambdaConL = -1*ones(dim,dim);

for i = 1:dim-1   
    if ((archiveL(i,:))~=-1)
        fit2 = archiveL(i,:);
    else
        x2 = x1;
        x2(i) = xL(i);
        Pop2 = Global.Evaluate(INDIVIDUAL({x2(1:u_dim),x2(u_dim+1:dim)}),type);
        fit2 =[Pop2.([type,'_obj']),Pop2.([type,'_con'])];
        archiveL(i,:) = fit2;
    end
    
    for j = u_dim+1:dim
        if j>i
            if ((archiveL(j,:))~=-1)
                fit3 = archiveL(j,:);
            else
                x3 = x1;
                x3(j) = xL(j);
                Pop3 = Global.Evaluate(INDIVIDUAL({x3(1:u_dim),x3(u_dim+1:dim)}),type);
                fit3 =[Pop3.([type,'_obj']),Pop3.([type,'_con'])];
                archiveL(j,:) = fit3;
            end
            x4 = x1;
            x4(i) = xL(i);
            x4(j) = xL(j);
            Pop4 = Global.Evaluate(INDIVIDUAL({x4(1:u_dim),x4(u_dim+1:dim)}),type);
            fit4 =[Pop4.([type,'_obj']),Pop4.([type,'_con'])];
            d1 = fit2(1:Global.M(2))-fit1(1:Global.M(2));
            d2 = fit4(1:Global.M(2))-fit3(1:Global.M(2));
            d3 = abs(fit2(1+Global.M(2):end)-fit1(1+Global.M(2):end));
            d4 = abs(fit4(1+Global.M(2):end)-fit2(1+Global.M(2):end));
            
            lambdaObjL(i,j,:) = abs(d1-d2);
            lambdaConL(i,j) = sum(d4.*d3)>0;
        end
    end
end

adjL = lambdaObjL(:,:,1)>epsilon | lambdaObjL(:,:,2)>epsilon | lambdaConL>0;

adj = adjL;
adj(logical(eye(dim))) = 1;
adj = adj|adj';
labels = findConnComp(adj);

groupNum = max(labels);
groups = cell(groupNum,1);

for i = 1:groupNum
    index = find(labels == i);
    if all(index <= u_dim) || all(index > u_dim)
        groups{i}{1} = index;
    else
        indexU = index(index <= u_dim);
        groups{i}{1} = indexU;
        indexL = index(index > u_dim);
        adjL1 = adjL(indexL,indexL);
        labels1 = findConnComp(adjL1);
        groupNum1 = max(labels1);
        for k = 1:groupNum1
            groups{i}{k+1} = [indexL(labels1==k)];
        end
    end
end
end

function labels = findConnComp(adj)

L = size(adj,1);

labels = zeros(1,L);
rts = [];
ccc = 0;

while true
    ind = find(labels==0);
    if ~isempty(ind)
        fue = ind(1);
        rts = [rts fue];
        list = [fue];
        ccc = ccc+1;
        labels(fue) = ccc;
        while true
            list_new = [];
            for lc = 1:length(list)
                p = list(lc);
                cp = find(adj(p,:));
                cx1 = cp(labels(cp)==0);
                labels(cx1)=ccc;
                list_new = [list_new cx1];
            end
            list = list_new;
            if isempty(list)
                break;
            end
        end
    else
        break;
    end
end
end
