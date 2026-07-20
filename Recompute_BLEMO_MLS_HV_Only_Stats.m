function Recompute_BLEMO_MLS_HV_Only_Stats()
% Recompute final upper-level HV statistics from saved BLEMO-MLS runs only.

    rootDir   = fileparts(mfilename('fullpath'));
    resultDir = fullfile(rootDir,'Results','BLEMO_MLS_FGTLEA_Comparison_7Problems_10Runs');
    runCsv    = fullfile(resultDir,'BLEMO_MLS_HV_only_recomputed_runs.csv');
    statsCsv  = fullfile(resultDir,'BLEMO_MLS_HV_only_median_iqr.csv');

    addpath(genpath(rootDir));

    problemRows = { ...
        'TP1',   'BLEMO_MLS_TP1_D1_2_R%d.mat'; ...
        'TP2',   'BLEMO_MLS_TP2_D1_14_R%d.mat'; ...
        'DS1_1', 'BLEMO_MLS_DS1_1_D10_10_R%d.mat'; ...
        'DS2_1', 'BLEMO_MLS_DS2_1_D10_10_R%d.mat'; ...
        'DS3_1', 'BLEMO_MLS_DS3_1_D10_10_R%d.mat'; ...
        'DS4',   'BLEMO_MLS_DS4_D1_9_R%d.mat'; ...
        'DS5',   'BLEMO_MLS_DS5_D1_9_R%d.mat'};

    runRows = struct([]);
    for p = 1:size(problemRows,1)
        problemName = problemRows{p,1};
        filePattern = problemRows{p,2};
        for r = 1:10
            matFile = fullfile(resultDir,sprintf(filePattern,r));
            row = emptyRunRow(problemName,r,matFile);
            try
                S = load(matFile,'Data');
                [row.HV,row.FeasibleCount,row.RangeCount,row.LegacyHV,row.LegacyCount] = ...
                    recomputeFinalHV(S.Data);
                row.Status = "OK";
            catch err
                row.Status = "ERROR";
                row.ErrorID = string(err.identifier);
                row.ErrorMessage = string(err.message);
            end
            runRows = [runRows; row]; %#ok<AGROW>
        end
    end

    RunHV = struct2table(runRows);
    writetable(RunHV,runCsv);

    Stats = rowsToHVStats(RunHV);
    writetable(Stats,statsCsv);

    save(fullfile(resultDir,'BLEMO_MLS_HV_only_median_iqr.mat'), ...
        'RunHV','Stats','resultDir','runCsv','statsCsv');

    fprintf('Saved HV run values to %s\n',runCsv);
    fprintf('Saved HV median/IQR to %s\n',statsCsv);
end

function row = emptyRunRow(problemName,run,matFile)
    row = struct( ...
        'Problem',string(problemName), ...
        'Run',run, ...
        'HV',NaN, ...
        'FeasibleCount',NaN, ...
        'RangeCount',NaN, ...
        'LegacyHV',NaN, ...
        'LegacyCount',NaN, ...
        'Status',"MISSING", ...
        'ErrorID',"", ...
        'ErrorMessage',"", ...
        'MatFile',string(matFile));
end

function [hv,feasibleCount,rangeCount,legacyHV,legacyCount] = recomputeFinalHV(Data)
    hv = 0;
    feasibleCount = 0;
    rangeCount = 0;
    legacyHV = 0;
    legacyCount = 0;

    if isempty(Data.result) || isempty(Data.result(end).Population) || isempty(Data.PF)
        return;
    end

    Population = Data.result(end).Population;
    feasible = Population.upper_feasibles & Population.lower_feasibles;
    Population = Population(feasible);
    feasibleCount = length(Population);
    if isempty(Population)
        return;
    end

    PopObj = Population.upper_objs;
    PopObj = filterOutsidePFRange(PopObj,Data.PF);
    rangeCount = size(PopObj,1);
    if ~isempty(PopObj)
        hv = HV(PopObj,Data.PF);
    end

    legacyObj = PopObj;
    if ismember(class(Data.problem),{'TP1','TP2'})
        legacyObj = filterBeyondNearestPF(legacyObj,Data.PF);
    end
    legacyCount = size(legacyObj,1);
    if ~isempty(legacyObj)
        legacyHV = HV(legacyObj,Data.PF);
    end
end

function PopObj = filterOutsidePFRange(PopObj,PF)
    if isempty(PopObj) || isempty(PF)
        return;
    end
    fmax = max(PF,[],1);
    PopObj(any(PopObj>repmat(fmax,size(PopObj,1),1),2),:) = [];
end

function PopObj = filterBeyondNearestPF(PopObj,PF)
    if isempty(PopObj) || isempty(PF)
        return;
    end
    [~,nearest] = min(pdist2(PopObj,PF),[],2);
    betterThanPF = all(PopObj < PF(nearest,:),2);
    PopObj(betterThanPF,:) = [];
end

function Stats = rowsToHVStats(RunHV)
    problemNames = unique(RunHV.Problem,'stable');
    Problems = strings(numel(problemNames),1);
    Runs = zeros(numel(problemNames),1);
    HVMedian = NaN(numel(problemNames),1);
    HVIQR = NaN(numel(problemNames),1);

    for i = 1:numel(problemNames)
        current = RunHV(RunHV.Problem == problemNames(i) & RunHV.Status == "OK",:);
        values = current.HV(~isnan(current.HV));
        Problems(i) = problemNames(i);
        Runs(i) = numel(values);
        if ~isempty(values)
            HVMedian(i) = median(values);
            HVIQR(i) = iqr(values);
        end
    end

    Stats = table(Problems,Runs,HVMedian,HVIQR, ...
        'VariableNames',{'Problem','Runs','HVMedian','HVIQR'});
end
