function Population = tentative_explore(Population,beta,theta)
    Global = GLOBAL.GetObj();
    Population = Global.Evaluate(Population,'upper');
    
    gen = 1;
    N = length(Population);
    Initial_Objs = Population.upper_objs;
    Distance = pdist2(Initial_Objs,Initial_Objs); 
    D = max(Distance,[],2);
    Z  = min(Population.upper_objs,[],1);
    
    ideal_points = {};
    nadir_points = {};
    Objs = {};
    
    MaxGen = floor(0.5*Global.maxFEs(1)/N);
    
    while gen <= MaxGen
%         
        P = [];
        for i=1:N
            P = cat(1,P,randperm(N,3));
        end
        Decs = DE_current_1_bin(Population(P(:,1)).decs,Population(P(:,2)).decs,Population(P(:,3)).decs,'bilevel',{1,0.5,0.1,20});
        Offspring = Global.Evaluate(INDIVIDUAL({Decs(:,1:Global.D(1)),Decs(:,1+Global.D(1):sum(Global.D))}),'upper');
        
%         Population = Select([Population,Offspring],N,W,'upper');
        Population = NSGAII_Update([Population,Offspring],N,'upper');
        
        labels = Population.labels(false);
        [Distance,ID] = min(pdist2(Population(labels).upper_objs,Initial_Objs),[],2);
        I = sum(Distance>D(ID))/length(Distance);
        
%         Global.upper_Output(Population);
        
        if I == 1
            break;
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