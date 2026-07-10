function [Score,PopObj] = HV(Population,optimum,problemName)
% <metric> <max>
% Hypervolume using the same final-population convention as FG-TLEA.

    if nargin < 3
        problemName = '';
    end

    if isa(optimum,'GLOBAL')
        Global = optimum;
        problemName = class(Global.problem);
        optimum = FGTLEAOptimumForHV(problemName,Global.PF);
    end

    [PopObj,ulCon,llCon] = ExtractPopulation(Population);
    if isempty(PopObj)
        Score = 0;
        return;
    end
    ulCon = NormalizeConstraints(ulCon,size(PopObj,1));
    llCon = NormalizeConstraints(llCon,size(PopObj,1));

    feasible = true(size(PopObj,1),1);
    if ~isempty(llCon)
        feasible = feasible & ~any(llCon > 0,2);
    end
    if ~isempty(ulCon)
        feasible = feasible & ~any(ulCon > 0,2);
    end
    PopObj = PopObj(feasible,:);

    if isempty(PopObj)
        Score = 0;
        return;
    end

    switch problemName
        case 'TP1'
            F2 = PopObj(:,2);
            F1 = -1 - F2 - sqrt(2*(F2+0.5).^2 + 0.5);
            PopObj(PopObj(:,1) < F1,:) = [];
        case 'TP2'
            F2 = PopObj(:,2);
            F1 = F2/2 + (1 - sqrt(F2/2)).^2;
            PopObj(PopObj(:,1) < F1,:) = [];
    end

    if isempty(PopObj)
        Score = 0;
        return;
    end

    [N,M] = size(PopObj);
    if size(optimum,2) ~= M
        Score = NaN;
        return;
    end
    fmin = min(min(PopObj,[],1),zeros(1,M));
    fmax = max(optimum,[],1);
    denominator = (fmax - fmin) * 1.1;

    PopObj = (PopObj - repmat(fmin,N,1))./repmat(denominator,N,1);
    PopObj(any(PopObj > 1,2),:) = [];

    RefPoint = ones(1,M) * 1.1;

    if isempty(PopObj)
        Score = 0;
    elseif M < 4
        pl = sortrows(PopObj);
        S  = {1,pl};
        for k = 1 : M-1
            S_ = {};
            for i = 1 : size(S,1)
                Stemp = Slice(cell2mat(S(i,2)),k,RefPoint);
                for j = 1 : size(Stemp,1)
                    temp(1) = {cell2mat(Stemp(j,1))*cell2mat(S(i,1))};
                    temp(2) = Stemp(j,2);
                    S_      = Add(temp,S_);
                end
            end
            S = S_;
        end
        Score = 0;
        for i = 1 : size(S,1)
            p     = Head(cell2mat(S(i,2)));
            Score = Score + cell2mat(S(i,1))*abs(p(M)-RefPoint(M));
        end
    else
        SampleNum = 1000000;
        MaxValue  = RefPoint;
        MinValue  = min(PopObj,[],1);
        Samples   = unifrnd(repmat(MinValue,SampleNum,1),repmat(MaxValue,SampleNum,1));
        for i = 1 : size(PopObj,1)
            drawnow();
            domi = true(size(Samples,1),1);
            m    = 1;
            while m <= M && any(domi)
                domi = domi & PopObj(i,m) <= Samples(:,m);
                m    = m + 1;
            end
            Samples(domi,:) = [];
        end
        Score = prod(MaxValue-MinValue)*(1-size(Samples,1)/SampleNum);
    end
end

function optimum = FGTLEAOptimumForHV(problemName,fallback)
% HV only uses max(optimum,[],1). These rows are exactly
% max(Problem.GetOptimum(10000),[],1) in the FG-TLEA problem classes.
    switch problemName
        case 'TP1'
            optimum = [-1,0];
        case 'TP2'
            optimum = [1,0.5];
        case 'DS1_1'
            optimum = [1.1,1.1];
        case 'DS2_1'
            optimum = [0.8089777225965074,0.0136728251065964];
        case 'DS3_1'
            optimum = [1.3999618794836017,0.9999408794349228];
        case 'DS4'
            optimum = [1,2];
        case 'DS5'
            optimum = [1,1.8];
        otherwise
            optimum = fallback;
    end
end

function [PopObj,ulCon,llCon] = ExtractPopulation(Population)
    ulCon = [];
    llCon = [];

    if isempty(Population)
        PopObj = [];
    elseif isa(Population,'INDIVIDUAL')
        PopObj = Population.upper_objs;
        ulCon = Population.upper_cons;
        llCon = Population.lower_cons;
    elseif isnumeric(Population)
        PopObj = Population;
    else
        PopObj = ReadMember(Population,'ulObjs');
        if isempty(PopObj)
            PopObj = ReadMember(Population,'upper_objs');
        end
        ulCon = ReadMember(Population,'ulCon');
        if isempty(ulCon)
            ulCon = ReadMember(Population,'upper_cons');
        end
        llCon = ReadMember(Population,'llCon');
        if isempty(llCon)
            llCon = ReadMember(Population,'lower_cons');
        end
    end
end

function value = ReadMember(obj,name)
    value = [];
    try
        value = obj.(name);
    catch
        try
            value = cat(1,obj.(name));
        catch
            value = [];
        end
    end
end

function Con = NormalizeConstraints(Con,N)
    if isempty(Con)
        return;
    end
    if size(Con,1) ~= N && size(Con,2) == N
        Con = Con';
    end
    if size(Con,1) ~= N && numel(Con) == N
        Con = reshape(Con,N,[]);
    end
end

function S = Slice(pl,k,RefPoint)
    p  = Head(pl);
    pl = Tail(pl);
    ql = [];
    S  = {};
    while ~isempty(pl)
        ql  = Insert(p,k+1,ql);
        p_  = Head(pl);
        cell_(1,1) = {abs(p(k)-p_(k))};
        cell_(1,2) = {ql};
        S   = Add(cell_,S);
        p   = p_;
        pl  = Tail(pl);
    end
    ql = Insert(p,k+1,ql);
    cell_(1,1) = {abs(p(k)-RefPoint(k))};
    cell_(1,2) = {ql};
    S  = Add(cell_,S);
end

function ql = Insert(p,k,pl)
    flag1 = 0;
    flag2 = 0;
    ql    = [];
    hp    = Head(pl);
    while ~isempty(pl) && hp(k) < p(k)
        ql = [ql;hp];
        pl = Tail(pl);
        hp = Head(pl);
    end
    ql = [ql;p];
    m  = length(p);
    while ~isempty(pl)
        q = Head(pl);
        for i = k : m
            if p(i) < q(i)
                flag1 = 1;
            else
                if p(i) > q(i)
                    flag2 = 1;
                end
            end
        end
        if ~(flag1 == 1 && flag2 == 0)
            ql = [ql;Head(pl)];
        end
        pl = Tail(pl);
    end
end

function p = Head(pl)
    if isempty(pl)
        p = [];
    else
        p = pl(1,:);
    end
end

function ql = Tail(pl)
    if size(pl,1) < 2
        ql = [];
    else
        ql = pl(2:end,:);
    end
end

function S_ = Add(cell_,S)
    n = size(S,1);
    m = 0;
    for k = 1 : n
        if isequal(cell_(1,2),S(k,2))
            S(k,1) = {cell2mat(S(k,1))+cell2mat(cell_(1,1))};
            m = 1;
            break;
        end
    end
    if m == 0
        S(n+1,:) = cell_(1,:);
    end
    S_ = S;
end
