%Platform for Bilevel Opimaization.
function StartRun
    %% reset
    close('all');
    clear;
    clc;
    restoredefaultpath;
    cd(fileparts(mfilename('fullpath')));
    addpath(genpath(cd));
    

    %% setting
    Problems = {
        'DS1_1' [10 10] [40,10,400];...
%         'DS2_1' [10 10] [40,10,400];...
%         'DS3_1' [10 10] [40,10,400];...
%         'DS4' [1 9] [20,40,200];...
%         'DS5' [1 9] [20,40,200];...
%         'TP1' [1 2] [10,10,60];...
%         'TP2' [1 14] [10,30,300];...
%         'DS1_2' [10 10] [40,10,400];...
%         'DS2_2' [10 10] [40,10,400];...
%         'DS3_2' [10 10] [40,10,400];...
               };

    Algorithms = {
        
         'BLEMO_MLS';...
%         'NBLEMO';...
%         'HBLEMO';...
%         'stMOBEA';...
%         'cG_BLEMO';...
               };

    Metrics = {
        'UHV'...
        'UIGD'...
        'LGD'...
               };
    
    maxFEs = {
            [2.5e4 2e6];...
%             [2.5e4 2e6];...
%             [2.5e4 2e6];...
%             [1e4 5e5];...
%             [1e4 5e5];...
%             [3e3 5e5];...
%             [5e3 5e5];...
%             [3.5e4 3.5e6];...
%             [3.5e4 3.5e6];...
%             [3.5e4 3.5e6];...
        };

   
 
    Runs = 1;
    Fig.ShowFig = 1;
    Fig.Figtype = 1;
    isparallel = 0;
    
%     Runs = 15;
%     Save = 5;
%     Fig.ShowFig = 0;
%     Fig.Figtype = 1;
%     isparallel = 1;

    % Start
    Task = [];
    for i=1:size(Problems,1)
        if ismember(Problems{i,1},{'DS4','DS5','TP1','TP2'})
            Save = 2;
        else
            Save = 5;
        end
        for j=1:length(Algorithms)
            for k=1:Runs
                Inputs = {'-problem',Problems(i,:),'-algorithm',Algorithms(j),'-maxFEs',maxFEs{i},'-save',Save,...
                          '-Fig',Fig,'-CalMetric',Metrics,'-run',k};
                Task = cat(1,Task,Inputs);
            end
        end
    end
    
    if ~isparallel
        for k=1:size(Task,1)
            Global = GLOBAL(Task{k,:});
            Global.Start();
        end
    else
        parpool('local',12)
        parfor k=1:size(Task,1)
            Global = GLOBAL(Task{k,:},'-isparallel',isparallel);
            Global.Start();
        end
        delete(gcp);
    end

%     if Runs>1
%         GLOBAL.Data_dealing(Problems,Algorithms,Runs,Metrics,'std');
%     end

end