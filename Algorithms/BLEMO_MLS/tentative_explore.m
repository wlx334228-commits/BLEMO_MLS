function Population = tentative_explore(Population,redundant_dec,beta,theta)
    Global = GLOBAL.GetObj();
    Population = Global.Evaluate(Population,'upper');
         
    gen = 1;
    N = length(Population);
    Initial_Pop = Population;
    
    Z  = min(Population.upper_objs,[],1);
    
    ideal_points = {};
    nadir_points = {};
    Objs = {};
    
    MaxGen = 200;%floor(0.25*Global.maxFEs(1)/N);
    
    [Population,~,FrontNo,CrowdDis] = NSGAII_Update(Population,N,'upper');
    
    while gen <= MaxGen
        
        
        Target = TournamentSelection(2,N,FrontNo,-CrowdDis);
        P = [];
        for i=1:N
            P = cat(1,P,randperm(N,2));
        end
        Target_decs = Population(Target).decs;
        Decs = DE_current_1_bin(Target_decs,Population(P(:,1)).decs,Population(P(:,2)).decs,'bilevel',{0.9,0.5,1,20});
        Decs(:,redundant_dec) = Target_decs(:,redundant_dec);
        Offspring = Global.Evaluate(INDIVIDUAL({Decs(:,1:Global.D(1)),Decs(:,1+Global.D(1):sum(Global.D))}),'upper');
        
        [Population,~,FrontNo,CrowdDis] = NSGAII_Update([Population,Offspring],N,'upper');
        
        Distance = max(Population.upper_decs)-min(Population.upper_decs);
        
        Global.upper_Output(Population,3);
        
        if sum(Distance<0.1)==Global.D(1)
            Population = Initial_Pop;
            gen = 1;
        end
        
        ideal_points = cat(2,ideal_points,{min([Z;Offspring.upper_objs],[],1)});
        nadir_points = cat(2,nadir_points,{max(Population.upper_objs,[],1)});
        Objs = cat(2,Objs,{Population.upper_objs});
        
        if length(ideal_points)==beta && length(nadir_points)==beta
            max_change = calc_maxchange(ideal_points,nadir_points);
            RefPoint = max(cat(1,Objs{:}),[],1);
            HV = cell2mat(cellfun(@(x) CalHV(x,RefPoint), Objs, 'UniformOutput',false));
            E_u = (max(HV)-min(HV))/(max(HV)+min(HV));
            if max_change <= theta || E_u <= theta
                break;
            end
            ideal_points(1) = [];
            nadir_points(1) = [];
            Objs(1) = [];
        end
        gen = gen +1;
    end
end

function max_change = calc_maxchange(ideal_points,nadir_points)
    ideal_points = cat(1,ideal_points{:});
    nadir_points = cat(1,nadir_points{:});
    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    rz = abs(ideal_points(end,:) - ideal_points(1,:)) ./ max(abs(ideal_points(1,:)),delta_value);
    nrz = abs(nadir_points(end,:) - nadir_points(1,:)) ./ max(abs(nadir_points(1,:)),delta_value);
    max_change = max([rz, nrz]);
end