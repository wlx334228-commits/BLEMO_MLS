function Run_BLEMO_MLS_7Problems_10Runs
% Run BLEMO-MLS on seven benchmark problems for 10 independent runs.
% The script records LGD, UHV, UIGD, upper/lower FEs, and running time.

    close all;
    clc;

    root = fileparts(mfilename('fullpath'));
    restoredefaultpath;
    cd(root);
    addpath(genpath(root));

    algorithm = 'BLEMO_MLS';
    metrics   = {'UHV','UIGD','LGD'};
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

    rows = repmat(emptyRow(),0,1);

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
                row.Status = sprintf('ERROR:%s',err.identifier);
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
                    row.UHV         = fieldOrNaN(S.Result,'UHV');
                    row.UIGD        = fieldOrNaN(S.Result,'UIGD');
                else
                    row.Status = 'MISSING_RESULT_STRUCT';
                end
            elseif strcmp(row.Status,'OK')
                row.Status = 'MISSING_RESULT_FILE';
            end

            rows(end+1,1) = row; %#ok<AGROW>
            writeSummary(rows,csvFile,matFile,problems,metrics);
        end
    end

    fprintf('\nFinished. Summary saved to:\n%s\n%s\n',csvFile,matFile);
end

function row = emptyRow()
    row = struct( ...
        'Problem','', ...
        'Run',NaN, ...
        'LGD',NaN, ...
        'UHV',NaN, ...
        'UIGD',NaN, ...
        'FEsu',NaN, ...
        'FEl',NaN, ...
        'RunningTime',NaN, ...
        'MaxFEsu',NaN, ...
        'MaxFEl',NaN, ...
        'Status','', ...
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

function writeSummary(rows,csvFile,matFile,problems,metrics)
    Summary = struct2table(rows);
    writetable(Summary,csvFile);
    save(matFile,'Summary','problems','metrics');
end
