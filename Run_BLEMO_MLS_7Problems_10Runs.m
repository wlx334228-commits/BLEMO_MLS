function Run_BLEMO_MLS_7Problems_10Runs
% Run BLEMO-MLS on seven benchmark problems for 10 independent runs.
% The script records LGD, FG-TLEA-compatible HV, UIGD, upper/lower FEs,
% and running time.

    close all;
    clc;

    root = fileparts(mfilename('fullpath'));
    restoredefaultpath;
    cd(root);
    addpath(genpath(root));

    algorithm = 'BLEMO_MLS';
    metrics   = {'HV','UIGD','LGD'};
    runs      = 10;

    problems = {
        'TP1',   [1 2],   [10 10  60], [3e3   5e5]; ...
        'TP2',   [1 14],  [10 30 300], [5e3   5e5]; ...
        'DS1_1', [10 10], [40 10 400], [2.5e4 2e6]; ...
        'DS2_1', [10 10], [40 10 400], [2.5e4 2e6]; ...
        'DS3_1', [10 10], [40 10 400], [2.5e4 2e6]; ...
        'DS4',   [1 9],   [20 40 200], [1e4   5e5]; ...
        'DS5',   [1 9],   [20 40 200], [1e4   5e5]; ...
    };

    fig.ShowFig = 0;
    fig.Figtype = 1;

    outputDir = fullfile(root,'Results','BLEMO_MLS_7Problems_10Runs');
    if ~exist(outputDir,'dir')
        mkdir(outputDir);
    end
    csvFile = fullfile(outputDir,'BLEMO_MLS_7Problems_10Runs_summary.csv');
    matFile = fullfile(outputDir,'BLEMO_MLS_7Problems_10Runs_summary.mat');
    errorLogFile = fullfile(outputDir,'BLEMO_MLS_7Problems_10Runs_errors.log');

    rows = repmat(emptyRow(),0,1);
    initializeSummaryFiles(csvFile,matFile,errorLogFile,rows,problems,metrics);

    for p = 1:size(problems,1)
        problemName = problems{p,1};
        problemD    = problems{p,2};
        problemN    = problems{p,3};
        maxFEs      = problems{p,4};
        problemSpec = {problemName,problemD,problemN};

        if ismember(problemName,{'DS4','DS5','TP1','TP2'})
            saveGap = 2;
        else
            saveGap = 5;
        end

        for run = 1:runs
            row = emptyRow();
            row.Problem       = problemName;
            row.Run           = run;
            row.MaxFEsu       = maxFEs(1);
            row.MaxFEl        = maxFEs(2);
            row.Status        = 'OK';
            row.MatFile       = resultFile(root,algorithm,problemName,problemD,run);

            fprintf('Running %s on %s, run %d/%d, maxFEs=[%g,%g]\n', ...
                algorithm, problemName, run, runs, maxFEs(1), maxFEs(2));

            try
                rng(run,'twister');
                inputs = {'-problem',problemSpec, ...
                          '-algorithm',{algorithm}, ...
                          '-maxFEs',maxFEs, ...
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

            if strcmp(row.Status,'OK') && exist(row.MatFile,'file')
                S = load(row.MatFile,'Result');
                if isfield(S,'Result')
                    row.RunningTime = fieldOrNaN(S.Result,'runtime');
                    row.FEsu        = fieldOrNaN(S.Result,'upper_FEs');
                    row.FEl         = fieldOrNaN(S.Result,'lower_FEs');
                    row.LGD         = fieldOrNaN(S.Result,'LGD');
                    row.HV          = fieldOrNaN(S.Result,'HV');
                    row.UIGD        = fieldOrNaN(S.Result,'UIGD');
                else
                    row.Status = 'MISSING_RESULT_STRUCT';
                end
            elseif strcmp(row.Status,'OK')
                row.Status = 'MISSING_RESULT_FILE';
            end

            rows(end+1,1) = row; %#ok<AGROW>
            appendRowToCsv(row,csvFile);
            saveSummaryMat(rows,matFile,problems,metrics);
            fprintf('Saved task result to %s\n',csvFile);
        end
    end

    fprintf('\nFinished. Summary saved to:\n%s\n%s\n',csvFile,matFile);
end

function initializeSummaryFiles(csvFile,matFile,errorLogFile,rows,problems,metrics)
    if exist(csvFile,'file')
        delete(csvFile);
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
    fprintf(fid,'Problem,Run,LGD,HV,UIGD,FEsu,FEl,RunningTime,MaxFEsu,MaxFEl,Status,ErrorID,ErrorMessage,MatFile\n');
    clear cleaner;

    saveSummaryMat(rows,matFile,problems,metrics);
end

function row = emptyRow()
    row = struct( ...
        'Problem','', ...
        'Run',NaN, ...
        'LGD',NaN, ...
        'HV',NaN, ...
        'UIGD',NaN, ...
        'FEsu',NaN, ...
        'FEl',NaN, ...
        'RunningTime',NaN, ...
        'MaxFEsu',NaN, ...
        'MaxFEl',NaN, ...
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

function value = fieldOrNaN(s,fieldName)
    if isfield(s,fieldName) && ~isempty(s.(fieldName))
        value = s.(fieldName);
    else
        value = NaN;
    end
end

function appendRowToCsv(row,csvFile)
    fid = fopen(csvFile,'a');
    assert(fid>0,'RunScript:FileOpenError','Cannot open summary CSV: %s',csvFile);
    cleaner = onCleanup(@()fclose(fid));
    fprintf(fid,'%s,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%s,%s,%s,%s\n', ...
        csvText(row.Problem), ...
        row.Run, ...
        row.LGD, ...
        row.HV, ...
        row.UIGD, ...
        row.FEsu, ...
        row.FEl, ...
        row.RunningTime, ...
        row.MaxFEsu, ...
        row.MaxFEl, ...
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
    save(matFile,'Summary','rows','problems','metrics');
end

function Summary = rowsToTable(rows)
    variableNames = {'Problem','Run','LGD','HV','UIGD','FEsu','FEl', ...
        'RunningTime','MaxFEsu','MaxFEl','Status','ErrorID','ErrorMessage','MatFile'};
    if isempty(rows)
        Summary = table('Size',[0 numel(variableNames)], ...
            'VariableTypes',{'string','double','double','double','double','double','double', ...
            'double','double','double','string','string','string','string'}, ...
            'VariableNames',variableNames);
        return;
    end

    Summary = table( ...
        string({rows.Problem}'), ...
        [rows.Run]', ...
        [rows.LGD]', ...
        [rows.HV]', ...
        [rows.UIGD]', ...
        [rows.FEsu]', ...
        [rows.FEl]', ...
        [rows.RunningTime]', ...
        [rows.MaxFEsu]', ...
        [rows.MaxFEl]', ...
        string({rows.Status}'), ...
        string({rows.ErrorID}'), ...
        string({rows.ErrorMessage}'), ...
        string({rows.MatFile}'), ...
        'VariableNames',variableNames);
end
