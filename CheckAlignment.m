% clear all;

dddir = getenv('DBSDATADIR');

fixfile = [dddir 'COT_CamAlign_Fixes.xlsx'];

% AllSubjStr{1} = 'VM05'

datadir = getenv('COTNNDATA');
load([datadir 'SS.mat']);

% SS = {'VM20'};
% SS = {'VM06'}
% SS = {'VM20'}
SS = {'PD09'}

ndim = 2;
% ndim = 3;

% AllSubjStr = {'VM25'};

for subi = 1:length(SS)
    SSi = SS{subi}

    % Where are the AO .mat files?
    AOmatdir = [dddir '/' SSi '/AO_Case_Data/matfiles'];
    
    % Where are the task info .xlsx files?
    Taskdir = [dddir '/' SSi '/' SSi '_Task'];
    
    % Where are the dlc video timestamp .txt files?
    Kintimedir = [dddir '/' SSi '/' SSi '_Compressed_Video/' SSi];
    
    % Where are the dlc csvs?
    if ndim==2
        kdir = dir([dddir SSi '\' SSi '-*-*-*-*']);
        Kindir = [dddir SSi '\' kdir(1).name '\videos\'];
    elseif ndim==3
        Kindir = [dddir '/' SSi '/3D_CSV/'];
    end
    
    % Where is the Key file?
    Keyfilename = [dddir '/' SSi '/' SSi '_Key'];
    
    % Read in the Key file
    Key = tdfread(Keyfilename);
    
    nfiles = length(Key.AODepth);
    disp(['Found ' num2str(nfiles) ' total files']);
    
    % Discover AO mat file naming format
    AOmatfiles = dir([AOmatdir '/*.mat']);
    AO1 = AOmatfiles(1).name;

    for fii = 1:nfiles
        fii

        Warning = [];
        
        % Read in the fixes .xlsx file
        % Load fixes file and extract relevant rows
        fixnotes = readcell(fixfile);
        fixheaders = fixnotes(1,:);
        fixtab = cell2table(fixnotes(2:end,:));
        fixtab.Properties.VariableNames = fixheaders;
        subfii = strcmp(fixtab{:,'Subject'},SSi) & fixtab{:,'Filenum'}==fii;
        
        figure(1); clf;
        fprintf(['Doing session ' num2str(fii) ' (video session ' num2str(Key.Vid(fii)) ')\n']);
        % Load in the AO mat data
        AOfilestub = ['D' num2str(Key.AODepth(fii), '%0.3f') 'F' num2str(Key.F(fii),'%0.4i')];
        AOfilestub = dir([AOmatdir '/*' AOfilestub '.mat']);
        AOfilestub = AOfilestub.name;
        M = load([AOmatdir '/' AOfilestub]);

        if isfield(M,'CDIG_IN_1_Up')
            TTL1_times = M.CDIG_IN_1_Down/(1000*M.CDIG_IN_1_KHz) + M.CDIG_IN_1_TimeBegin;
        else
            TTL1_times = [];
        end
        if isfield(M,'CDIG_IN_2_Up')
            TTL2_times = M.CDIG_IN_2_Up/(1000*M.CDIG_IN_2_KHz) + M.CDIG_IN_2_TimeBegin;
        else
            TTL2_times = [];
        end
        hasthree = false;
        if isfield(M,'CDIG_IN_3_Up')
            hasthree = true;
            disp('3 TTL streams detected!');
            TTL3_times = M.CDIG_IN_3_Up/(1000*M.CDIG_IN_3_KHz) + M.CDIG_IN_3_TimeBegin;
        end
        
        % Make sure we know which TTL is which, store camera frame times as the
        % TTL receive times
        if ~hasthree
            if numel(TTL1_times) > numel(TTL2_times)
                camTime = M.CDIG_IN_1_Up/(1000*M.CDIG_IN_1_KHz) + M.CDIG_IN_1_TimeBegin;
                AO_Task_TTL_times = M.CDIG_IN_2_Down/(1000*M.CDIG_IN_2_KHz) + M.CDIG_IN_2_TimeBegin;
            else
                camTime = TTL2_times;
                AO_Task_TTL_times = TTL1_times;
            end
        else % THIS ASSUMES THAT IF THERE'S THREE TTLS, OUR CAMERA IS RECORDING AT A HIGHER FRAMERATE
            [~,mini] = min([median(diff(TTL1_times)), median(diff(TTL2_times)), median(diff(TTL3_times))]);
            [~,maxi] = max([median(diff(TTL1_times)), median(diff(TTL2_times)), median(diff(TTL3_times))]);
            eval(['camTime =  M.CDIG_IN_' num2str(mini) '_Up/(1000*M.CDIG_IN_' num2str(mini) '_KHz) + M.CDIG_IN_' num2str(mini) '_TimeBegin;']);
            eval(['AO_Task_TTL_times = M.CDIG_IN_' num2str(maxi) '_Down/(1000*M.CDIG_IN_' num2str(maxi) '_KHz) + M.CDIG_IN_' num2str(maxi) '_TimeBegin;']);
        end
        
        % Load in the task timing data. This is setup for VM01 - xxxx, add
        % something else for different format!
        X = readcell([Taskdir '/Session ' num2str(Key.Tas(fii)) '/Times_' SSi '-T' num2str(Key.Tas(fii)) '.xlsx']);

        % Clean it up
        X(1,:) = [];
        times = cell2mat(X(:,1));
        validtimes = find(times~=0);
        nvalid = length(validtimes);
        times = times(validtimes);
        namesdirty = X(validtimes, 3);
        
        namesclean = cell(nvalid,1);
        nameidx = nan(nvalid,1);
        for j = 1:length(validtimes)
            splitname = strsplit(namesdirty{j},' ');
            namesclean{j} = splitname{1};
            if length(splitname) > 1
                nameidx(j) = str2double(splitname{end});
            end
        end
        
        xind = X(validtimes, 4);
        for xindi = 1:length(xind)
            if ismissing(xind{xindi})
                xind{xindi} = NaN;
            end
        end
        nameidx = cell2mat(xind);
        
        % Do the TTL regression to convert task timestamps to AO times
        if any(strcmp(SSi, {'PD03', 'PD02', 'PD05', 'VM03'}))
            if any(strcmp(namesclean, 'CameraFrame'))
                Task_TTL_times = times(strcmp(namesclean, 'CameraFrame'));
            else
                Task_TTL_times = times(strcmp(namesclean, 'BeforeStimuliShown'));
            end
        else
            Task_TTL_times = times(strcmp(namesclean,'TTL'));
        end
% %         
%         if fii == 2
%             plot(diff(AO_Task_TTL_times), 'r', 'Marker', '.');
%             hold on;
%             plot(diff(Task_TTL_times), 'b');
%             STOPP
%         end

        % DO AUTO FIXING OF TASK TTLS
        if ~isempty(subfii)
            if ~isnan(fixtab{subfii,'Task_TTL_Start'})
                disp('Applied Task TTL trim');
                if isnan(fixtab{subfii,'Task_TTL_Stop'})
                    Task_TTL_times = Task_TTL_times(fixtab{subfii,'Task_TTL_Start'}:end);
                elseif sign(fixtab{subfii,'Task_TTL_Stop'}) == -1
                    Task_TTL_times = Task_TTL_times(fixtab{subfii,'Task_TTL_Start'}:end+fixtab{subfii,'Task_TTL_Stop'});
                else
                    Task_TTL_times = Task_TTL_times(fixtab{subfii,'Task_TTL_Start'}:fixtab{subfii,'Task_TTL_Stop'});
                end
            end
            if ~isnan(fixtab{subfii,'AO_TTL_Start'})
                disp('Applied AO TTL trim');
                AO_Task_TTL_times = AO_Task_TTL_times(fixtab{subfii,'AO_TTL_Start'}:end);
            end
        end
        origAO = AO_Task_TTL_times;
        % Try to auto-clip bad AO task TTLs out
        % First clip the end if it's bad
        badAO = find(diff(origAO) > 15, 1, 'last');
        if badAO == length(origAO)-1
            origAO = origAO(1:end-1);
            disp('Clipped last TTL');
            Warning = [Warning ' Clipped last AO TTL '];
        end
        % Then clip the start if it's bad
        badAO = find(diff(origAO) > 15, 1, 'last');
        if badAO < length(origAO)/2
            origAO = origAO(badAO+1:end);
            Warning = [Warning ' Clipped Start AO TTL '];
        end
        newAO = origAO;
        
        % If there are still more received than sent task TTLs, try to scan
        % for them
        ldif = length(newAO) - length(Task_TTL_times);
        if ldif > 0
            Warning = [Warning ' Did regression TTL match '];
            ntry = ldif+1;
            trycorr = nan(ntry,1);
            for tryi = 1:ntry
                trycorr(tryi) = corr(diff(newAO((1:length(Task_TTL_times))+tryi-1))', diff(Task_TTL_times));
            end
            [~,maxi] = max(trycorr);
            newAO = newAO((1:length(Task_TTL_times))+maxi-1);
        else
            Warning = [Warning ' TTL nums matched '];
        end
        AO_Task_TTL_times = newAO;
        
        [b, ~, ~, ~, stats] = regress(AO_Task_TTL_times', [Task_TTL_times, ones(length(Task_TTL_times),1)]);
        if stats(1) < 0.999999
            disp('Potential problem with Task TTLs, CHECK PLOTS! Moving to next session.');
            disp(stats(1))
            disp(b)
            plot(diff(AO_Task_TTL_times), 'r');
            hold on;
            plot(diff(Task_TTL_times), 'b');
            break;
        end

        % Plot
        figure(1);
%         figure;
        subplot(2,2,1);
        TaskConv = Task_TTL_times*b(1) + b(2);
        plot(AO_Task_TTL_times, 'r', 'Marker', '.');
        hold on;
        plot(TaskConv,'b');
        
        subplot(2,2,2);
        plot(diff(AO_Task_TTL_times), 'r', 'Marker', '.');
        hold on;
        plot(diff(Task_TTL_times), 'b');
        
        % Rectify timestamps for kinematics...
        topvid = dir([Kintimedir '/session' num2str(Key.Vid(fii),'%0.3i') '/*_' SSi '*topCam*.mp4']);
        datestr = topvid(1).name(1:8);
        
        kinfiles = dir([Kindir '*.csv']);
        kinsplit = strsplit(kinfiles(1).name,'session');
        DLC = readcell([Kindir '/' kinsplit{1} 'session' num2str(Key.Vid(fii),'%0.3i') kinsplit{2}(4:end)]);
        if ndim == 2
            nframe = max(cell2mat(DLC(4:end,1)))+1;
        elseif ndim == 3
            nframe = max(cell2mat(DLC(2:end,end)))+1;
        end
        
        kintimes = cell2mat(readcell([Kintimedir '/session' num2str(Key.Vid(fii),'%0.3i') '/' datestr '_' SSi '_session' num2str(Key.Vid(fii),'%0.3i') '_topCam_timestamps.txt']));
        kintimes = kintimes(1:2:end)/(1e9);

        pauseAlign = 0;
        didFix = 0;
        
        camAligni = find(diff(camTime)>1,1,'first');
        kinAligni = find(kintimes>1,1,'first');

        % if strcmp(SSi, 'PD02') && fii == 2
        %     STOPP
        % end
        % Apply fixes to alignment indices
        if ~isempty(subfii)
            if ~isnan(fixtab{subfii,'pauseAlign'})
                pauseAlign = fixtab{subfii,'pauseAlign'};
                disp('Doing pause alignment');
            else
                if ~isnan(fixtab{subfii,'camAligni_direct'})
                    camAligni = fixtab{subfii,'camAligni_direct'};
                    disp('Did cam align index fix');
                    didFix = 1;
                end
                if ~isnan(fixtab{subfii,'camAligni_gap'})
                    gaplocs = find(diff(camTime)>1);
                    camAligni = gaplocs(fixtab{subfii,'camAligni_gap'}) + fixtab{subfii,'camAligni_gap_adj'};
                    disp('Did cam align index fix');
                    didFix = 1;
                end
                if ~isnan(fixtab{subfii,'kinAligni'})
                    kinAligni = fixtab{subfii,'kinAligni'};
                    disp('Did kin align index fix');
                    didFix = 1;
                end
            end
        end

        if pauseAlign
            camAligni = find(diff(camTime)>0.95 & diff(camTime)<1.05);
            kinAligni = find(kintimes>0.95 & kintimes<1.05);
        end

        % FIX FOR BAD REPEAT FIRST FRAME!!!
        % Adjust so that we are matching the samples AFTER the delay
        camAligni = camAligni+1;
        kintimes = kintimes(2:end);

        ktadd = nan(nframe,1);
        ktadd(1:kinAligni-1) = flipud(-cumsum(flipud(kintimes(1:kinAligni-1))));
        ktadd(kinAligni) = 0;
        % SPECIAL FIX BECAUSE IT APPEARS THAT THE LAST FEW SAMPLES GET CUT
        % OFF IN THE 3D CSV FILES???
%         ktadd(kinAligni+1:end) = cumsum(kintimes(kinAligni:end-1));
        ktadd(kinAligni+1:end) = cumsum(kintimes(kinAligni:end-(length(kintimes)-length(ktadd)+1)));
        K.Time = camTime(camAligni)+ktadd;

        % How many frames before aligni?
        if camAligni - (kinAligni-1) > 0
            camClip = camTime(camAligni-(kinAligni-1):end);
        else
            camClip = [nan(1,kinAligni-camAligni) camTime];
        end

        plotKtime = K.Time;
               
        figure(1);
        subplot(2,2,3);
        plot(camClip, 'r', 'Marker', '.');
        hold on;
        plot(plotKtime, 'b')
        
        figure(1);
        subplot(2,2,4);
        plot(diff(camClip), 'r', 'Marker', '.');
        hold on;
        plot(diff(plotKtime), 'b');
        
        if didFix
            Warning = [Warning ' DID FIX'];
        end
        if pauseAlign
            Warning = [Warning ' PAUSE ALIGNED'];
        end
        if nframe > numel(camTime)
            Warning = [Warning ' CAM WARNING'];
        end
        sgtitle([SSi ' ' num2str(fii) '/' num2str(nfiles) Warning]);

        camAligni
        kinAligni

        while ~waitforbuttonpress
        end
        
%         plot(K.Time,1:length(K.Time), '.r')

    end
end
close all;