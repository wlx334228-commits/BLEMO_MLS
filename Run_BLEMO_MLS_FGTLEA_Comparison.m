function Run_BLEMO_MLS_FGTLEA_Comparison
% Run BLEMO-MLS under the comparison settings used for FG-TLEA.
% The compared algorithm population size is Nu*Nl for both levels.

    close all;
    clc;

    root = fileparts(mfilename('fullpath'));
    restoredefaultpath;
    cd(root);
    addpath(genpath(root));

    algorithm = 'BLEMO_MLS';
    metrics   = {'UIGD','LGD'};
    runs      = 10;

    problems = {
        'TP1',   [1 2],   5,  12,  60, 5e5; ...
        'TP2',   [1 14],  5,  60, 300, 5e5; ...
        'DS1_1', [10 10], 20, 20, 400, 1.5e6; ...
        'DS2_1', [10 10], 20, 20, 400, 3e6; ...
        'DS3_1', [10 10], 20, 20, 400, 3e6; ...
        'DS4',   [1 9],   5,  40, 200, 1.5e6; ...
        'DS5',   [1 9],   5,  40, 200, 1.5e6; ...
    };

    fig.ShowFig = 0;
    fig.Figtype = 1;

    outputDir = fullfile(root,'Results','BLEMO_MLS_FGTLEA_Comparison_7Problems_10Runs');
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    csvFile      = fullfile(outputDir,'BLEMO_MLS_FGTLEA_Comparison_summary.csv');
    statsCsvFile = fullfile(outputDir,'BLEMO_MLS_FGTLEA_Comparison_stats.csv');
    matFile      = fullfile(outputDir,'BLEMO_MLS_FGTLEA_Comparison_summary.mat');
    errorLogFile = fullfile(outputDir,'BLEMO_MLS_FGTLEA_Comparison_errors.log');

    rows = repmat(emptyRow(),0,1);
    initializeSummaryFiles(csvFile,statsCsvFile,matFile,errorLogFile,rows,problems,metrics);

    for p = 1:size(problems,1)
        problemName = problems{p,1};
        problemD    = problems{p,2};
        Nu          = problems{p,3};
        Nl          = problems{p,4};
        upperLowerN = Nu*Nl;
        eliteN      = problems{p,5};
        problemN    = [upperLowerN upperLowerN eliteN];
        maxFE       = problems{p,6};
        problemSpec = {problemName,problemD,problemN};

        if ismember(problemName,{'DS4','DS5','TP1','TP2'})
            saveGap = 2;
        else
            saveGap = 5;
        end

        for run = 1:runs
            row = emptyRow();
            row.Problem = problemName;
            row.Run     = run;
            row.Nu      = Nu;
            row.Nl      = Nl;
            row.NUpperLower = upperLowerN;
            row.NElite      = eliteN;
            row.MaxFE   = maxFE;
            row.Status  = 'OK';

            frameworkMatFile = resultFile(root,algorithm,problemName,problemD,run);
            row.MatFile = archivedResultFile(outputDir,algorithm,problemName,problemD,run);

            fprintf('Running %s on %s, run %d/%d, Nu=%d, Nl=%d, N=[%d,%d,%d], maxFE=%g\n', ...
                algorithm, problemName, run, runs, Nu, Nl, problemN(1), problemN(2), problemN(3), maxFE);

            runStartedAt = datetime('now');
            try
                rng(run,'twister');
                inputs = {'-problem',problemSpec, ...
                          '-algorithm',{algorithm}, ...
                          '-maxFEs',maxFE, ...
                          '-save',saveGap, ...
                          '-Fig',fig, ...
                          '-CalMetric',metrics, ...
                          '-run',run};
                Global = GLOBAL(inputs{:});
                Global.Start();
            catch err
                row.Status = 'ERROR';
                row.ErrorID = err.identifier;
                row.ErrorMessage = oneLine(err.message);
                appendErrorLog(err,errorLogFile,problemName,run);
                warning('BLEMO-MLS run failed on %s R%d: %s', ...
                    problemName,run,err.message);
            end

            if strcmp(row.Status,'OK') && exist(frameworkMatFile,'file') && resultWasUpdated(frameworkMatFile,runStartedAt)
                S = load(frameworkMatFile,'Data','Result');
                if isfield(S,'Result')
                    row.RunningTime = fieldOrNaN(S.Result,'runtime');
                    row.FEsu        = fieldOrNaN(S.Result,'upper_FEs');
                    row.FEl         = fieldOrNaN(S.Result,'lower_FEs');
                    row.TotalFE     = row.FEsu + row.FEl;
                    row.LGD         = fieldOrNaN(S.Result,'LGD');
                    if isfield(S,'Data')
                        row.HV = finalFGTLEAHV(S.Data);
                    end
                    row.UIGD        = fieldOrNaN(S.Result,'UIGD');
                    copyfile(frameworkMatFile,row.MatFile);
                else
                    row.Status = 'MISSING_RESULT_STRUCT';
                end
            elseif strcmp(row.Status,'OK')
                if exist(frameworkMatFile,'file')
                    row.Status = 'STALE_RESULT_FILE';
                    row.ErrorMessage = sprintf('Result file was not updated by this run: %s',frameworkMatFile);
                else
                    row.Status = 'MISSING_RESULT_FILE';
                    row.ErrorMessage = sprintf('Missing result file: %s',frameworkMatFile);
                end
            end

            rows(end+1,1) = row; %#ok<AGROW>
            appendRowToCsv(row,csvFile);
            saveSummaryMat(rows,matFile,problems,metrics);
            writeStatsCsv(rows,statsCsvFile);
            fprintf('Saved summary row to %s\n',csvFile);
        end
    end

    fprintf('\nFinished. Summary saved to:\n%s\n%s\n',csvFile,matFile);
end

function initializeSummaryFiles(csvFile,statsCsvFile,matFile,errorLogFile,rows,problems,metrics)
    if exist(csvFile,'file')
        delete(csvFile);
    end
    if exist(statsCsvFile,'file')
        delete(statsCsvFile);
    end
    if exist(matFile,'file')
        delete(matFile);
    end
    if exist(errorLogFile,'file')
        delete(errorLogFile);
    end

    fid = fopen(csvFile,'w');
    assert(fid>0,'RunScript:FileOpenError','Cannot open summary CSV: %s',csvFile);
    cleaner = onCleanup(@()fclose(fid));
    fprintf(fid,'Problem,Run,Nu,Nl,NUpperLower,NElite,LGD,HV,UIGD,FEsu,FEl,TotalFE,RunningTime,MaxFE,Status,ErrorID,ErrorMessage,MatFile\n');
    clear cleaner;

    saveSummaryMat(rows,matFile,problems,metrics);
    writeStatsCsv(rows,statsCsvFile);
end

function row = emptyRow()
    row = struct( ...
        'Problem','', ...
        'Run',NaN, ...
        'Nu',NaN, ...
        'Nl',NaN, ...
        'NUpperLower',NaN, ...
        'NElite',NaN, ...
        'LGD',NaN, ...
        'HV',NaN, ...
        'UIGD',NaN, ...
        'FEsu',NaN, ...
        'FEl',NaN, ...
        'TotalFE',NaN, ...
        'RunningTime',NaN, ...
        'MaxFE',NaN, ...
        'Status','', ...
        'ErrorID','', ...
        'ErrorMessage','', ...
        'MatFile','');
end

function file = resultFile(root,algorithm,problemName,problemD,run)
    file = fullfile(root,'Data',algorithm, ...
        sprintf('%s_%s_D%d_%d_R%d.mat',algorithm,problemName, ...
        problemD(1),problemD(2),run));
end

function file = archivedResultFile(outputDir,algorithm,problemName,problemD,run)
    file = fullfile(outputDir, ...
        sprintf('%s_%s_D%d_%d_R%d.mat',algorithm,problemName, ...
        problemD(1),problemD(2),run));
end

function value = fieldOrNaN(s,fieldName)
    if isfield(s,fieldName) && ~isempty(s.(fieldName))
        value = s.(fieldName);
    else
        value = NaN;
    end
end

function score = finalFGTLEAHV(Data)
    score = NaN;
    if isempty(Data.result) || isempty(Data.result(end).Population)
        return;
    end

    PopObj = finalFeasibleUpperObjs(Data.result(end).Population,Data.PF);
    if isempty(PopObj) || isempty(Data.PF)
        score = 0;
        return;
    end

    score = HV(PopObj,Data.PF);
end

function PopObj = finalFeasibleUpperObjs(Population,PF)
    PopObj = [];
    if isempty(Population)
        return;
    end

    feasible = Population.upper_feasibles & Population.lower_feasibles;
    Population = Population(feasible);
    if isempty(Population)
        return;
    end

    PopObj = Population.upper_objs;
    PopObj = filterOutsidePFRange(PopObj,PF);
end

function PopObj = filterOutsidePFRange(PopObj,PF)
    if isempty(PopObj) || isempty(PF)
        return;
    end
    fmax = max(PF,[],1);
    PopObj(any(PopObj>repmat(fmax,size(PopObj,1),1),2),:) = [];
end

function tf = resultWasUpdated(file,runStartedAt)
    info = dir(file);
    tf = ~isempty(info) && datetime(info.datenum,'ConvertFrom','datenum') >= runStartedAt;
end

function appendRowToCsv(row,csvFile)
    fid = fopen(csvFile,'a');
    assert(fid>0,'RunScript:FileOpenError','Cannot open summary CSV: %s',csvFile);
    cleaner = onCleanup(@()fclose(fid));
    fprintf(fid,'%s,%d,%d,%d,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%s,%s,%s,%s\n', ...
        csvText(row.Problem), ...
        row.Run, ...
        row.Nu, ...
        row.Nl, ...
        row.NUpperLower, ...
        row.NElite, ...
        row.LGD, ...
        row.HV, ...
        row.UIGD, ...
        row.FEsu, ...
        row.FEl, ...
        row.TotalFE, ...
        row.RunningTime, ...
        row.MaxFE, ...
        csvText(row.Status), ...
        csvText(row.ErrorID), ...
        csvText(row.ErrorMessage), ...
        csvText(row.MatFile));
    clear cleaner;
end

function appendErrorLog(err,errorLogFile,problemName,run)
    fid = fopen(errorLogFile,'a');
    assert(fid>0,'RunScript:FileOpenError','Cannot open error log: %s',errorLogFile);
    cleaner = onCleanup(@()fclose(fid));
    fprintf(fid,'\n==== %s run %d ====\n',problemName,run);
    fprintf(fid,'%s\n',getReport(err,'extended','hyperlinks','off'));
    clear cleaner;
end

function value = oneLine(value)
    value = strtrim(regexprep(char(value),'\s+',' '));
end

function value = csvText(value)
    value = char(value);
    value = strrep(value,'"','""');
    value = ['"',value,'"'];
end

function saveSummaryMat(rows,matFile,problems,metrics)
    Summary = rowsToTable(rows);
    Stats = rowsToStats(rows);
    save(matFile,'Summary','Stats','rows','problems','metrics');
end

function writeStatsCsv(rows,statsCsvFile)
    Stats = rowsToStats(rows);
    writetable(Stats,statsCsvFile);
end

function Summary = rowsToTable(rows)
    variableNames = {'Problem','Run','Nu','Nl','NUpperLower','NElite','LGD','HV','UIGD', ...
        'FEsu','FEl','TotalFE','RunningTime','MaxFE','Status','ErrorID', ...
        'ErrorMessage','MatFile'};
    if isempty(rows)
        Summary = table('Size',[0 numel(variableNames)], ...
            'VariableTypes',{'string','double','double','double','double','double', ...
            'double','double','double','double','double','double','double', ...
            'double','string','string','string','string'}, ...
            'VariableNames',variableNames);
        return;
    end

    Summary = table( ...
        string({rows.Problem}'), ...
        [rows.Run]', ...
        [rows.Nu]', ...
        [rows.Nl]', ...
        [rows.NUpperLower]', ...
        [rows.NElite]', ...
        [rows.LGD]', ...
        [rows.HV]', ...
        [rows.UIGD]', ...
        [rows.FEsu]', ...
        [rows.FEl]', ...
        [rows.TotalFE]', ...
        [rows.RunningTime]', ...
        [rows.MaxFE]', ...
        string({rows.Status}'), ...
        string({rows.ErrorID}'), ...
        string({rows.ErrorMessage}'), ...
        string({rows.MatFile}'), ...
        'VariableNames',variableNames);
end

function Stats = rowsToStats(rows)
    variableNames = {'Problem','Runs','HVMedian','HVIQR', ...
        'RunningTimeMedian','RunningTimeIQR'};
    if isempty(rows)
        Stats = table('Size',[0 numel(variableNames)], ...
            'VariableTypes',{'string','double','double','double','double','double'}, ...
            'VariableNames',variableNames);
        return;
    end

    problemNames = unique({rows.Problem},'stable');
    Problems = strings(numel(problemNames),1);
    Runs = zeros(numel(problemNames),1);
    HVMedian = NaN(numel(problemNames),1);
    HVIQR = NaN(numel(problemNames),1);
    RunningTimeMedian = NaN(numel(problemNames),1);
    RunningTimeIQR = NaN(numel(problemNames),1);

    for i = 1:numel(problemNames)
        problem = problemNames{i};
        inProblem = strcmp({rows.Problem},problem);
        okRows = inProblem & strcmp({rows.Status},'OK');
        current = rows(okRows);

        Problems(i) = string(problem);
        Runs(i) = numel(current);
        HVMedian(i) = medianOrNaN([current.HV]);
        HVIQR(i) = iqrOrNaN([current.HV]);
        RunningTimeMedian(i) = medianOrNaN([current.RunningTime]);
        RunningTimeIQR(i) = iqrOrNaN([current.RunningTime]);
    end

    Stats = table(Problems,Runs,HVMedian,HVIQR,RunningTimeMedian, ...
        RunningTimeIQR,'VariableNames',variableNames);
end

function value = medianOrNaN(values)
    values = values(~isnan(values));
    if isempty(values)
        value = NaN;
    else
        value = median(values);
    end
end

function value = iqrOrNaN(values)
    values = values(~isnan(values));
    if isempty(values)
        value = NaN;
    else
        value = iqr(values);
    end
end
