% Create data tartare: Raw-ish data presented in a tasty form
%
% Version 2
%
% Requires the following pre-processed data files:
% 1. Converted .mat files from the AlphaOmega .mpx files
% 2. Sorted .plx files identifying the spikes from the AO data
% 3. Task output .csv files identifying the task timestamps and targets
% 4. DLC kinematic .csv files resulting from DLC-analyzed videos
% 5. A "Key" file (tab-separated) which denotes which sessions go together
%
% Requires that the Plexon Matlab SDK is on your Matlab path
%A
% This is a preliminary version that works for early subjects (PD01-03,
% VM01). It's a bit cumbersome due to changing data format standards.
%
% This version expects 3D kinematics from DLC analysis of a single camera.
%
% Outputs tasty .mat files with spike times, task times and kinematics, all
% synchronized into the AO clock. One file is created for each session, and
% can contain multiple spikes from multiple channels.
%
% Different data are in different structures.
% N: Contains the neural data (spike times, spike IDs)
% T: Contains the task data
% K: Contains the kinematic data
%
% All data fields are explained in detail in Data_Tartare_v0_notes.xlsx
%
% Rex Tien, July 2021
%
% Updated Sept 2022 to utilize 3D kinematic data
% Updated May 2023 to save files with no sorted units
% Updated Jan 2024 to be a function and to utilize a settings file for
%   camera alignment fixes, instead of hardcoding them here in the script.
%   Also saving the LFP and SPK fields from the AO file

function create_data_tartare(dddir, SSi, savedir, fixfile, ndim)

disp(['Doing ' SSi ]);

% Where are the AO .mat files?
AOmatdir = [dddir '/' SSi '/AO_Case_Data/matfiles'];

% Where are the sorted AO .plx files?
AOplxdir = ['C:/Users/RexBoxOne/Desktop/AOPlex/' SSi '/Sorted'];

% Where are the task info .xlsx files?
Taskdir = [dddir '/' SSi '/' SSi '_Task'];

% Where are the kinematic .csv files?
% Kindir = ['C:/Users/RexBoxOne/Documents/DeepLabCut/dbsdlc/' SSi '/COTask-Rex-2021-07-07/videos'];
if ndim == 2
    Kindir = dir([dddir '/' SSi '/' SSi '-*-*-*-*']);
    Kindir = [dddir '/' SSi '/' Kindir(1).name '/videos/'];
elseif ndim == 3
    Kindir = [dddir '/' SSi '/3D_CSV/'];
end

% Where are the dlc video timestamp .txt files?
Kintimedir = [dddir '/' SSi '/' SSi '_Compressed_Video/' SSi];

% Where is the Key file?
Keyfilename = [dddir '/' SSi '/' SSi '_Key'];

% Read in the Key file
Key = tdfread(Keyfilename);

nfiles = length(Key.AODepth);

for fii = 1:nfiles
    disp(['Doing file ' num2str(fii) ' (Depth ' num2str(Key.AODepth(fii)) ')']);
    
    % Read in the fixes .xlsx file
    % Load fixes file and extract relevant rows
    fixnotes = readcell(fixfile);
    fixheaders = fixnotes(1,:);
    fixtab = cell2table(fixnotes(2:end,:));
    fixtab.Properties.VariableNames = fixheaders;
    subfii = strcmp(fixtab{:,'Subject'},SSi) & fixtab{:,'Filenum'}==fii;
    
    % Load in the AO mat data
    AOfileload = dir([AOmatdir '/*D' num2str(Key.AODepth(fii), '%0.3f') 'F' num2str(Key.F(fii),'%0.4i') '*.mat']);
    M = load([AOmatdir '/' AOfileload.name]);
    
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
    
    % Store the depth of the recording
    N.Depth = Key.AODepth(fii);
    
    % Now to get the sorted spike data.
    % First determine how many chanels we have with sorted units
    if exist(AOplxdir ,'dir')
        sortfi = dir([AOplxdir '/*D' num2str(Key.AODepth(fii), '%0.3f') 'F' num2str(Key.F(fii),'%0.4i') '*.plx']);
        sortfi = {sortfi.name};
        if isempty(sortfi)
            disp('No units found - creating empty N file.');
        end
        nsortfi = length(sortfi);
        N.SpkID = [];
        itotspike = 1;
        for isortfi = 1:nsortfi
            chanstr = extractBetween(sortfi{isortfi},'_','.plx');
            chanstr = chanstr{1};
            channo = round(str2double(chanstr));
        
            % Read the plx file
            [OpenedFileName, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = plx_information([AOplxdir '/' sortfi{isortfi}]);
            % Read units until there are no more, excluding unsorted
            for uni = 1:10
                [nts, ts] = plx_ts(OpenedFileName, channo , uni );
                if nts == 0
                    break;
                end
                N.SpkTimes{itotspike} = ts + M.(['CSPK_' num2str(channo,'%0.2i') '_TimeBegin']);
                N.SpkID = [N.SpkID; [channo, uni]];
                itotspike = itotspike+1;
            end
        end
        if isempty(sortfi)
            N.SpkTimes = cell(0);
        end
    end
    
    % Now also save the CLFP and CSPK data in N
    fnm = fieldnames(M);
    clfpnames = fnm(contains(fnm,'CLFP'));
    cspknames = fnm(contains(fnm,'CSPK'));
    for lfpi = 1:length(clfpnames)
        N.LFP.(clfpnames{lfpi}) = M.(clfpnames{lfpi});
    end
    for spki = 1:length(cspknames)
        N.SPK.(cspknames{spki}) = M.(cspknames{spki});
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
    
    % Do the TTL regression to convert task timestamps to AO times. VM03 didn't
    % have TTL named timestamp, so use CameraFrame
    if any(strcmp(SSi, {'PD03', 'PD02', 'PD05', 'VM03'}))
        if any(strcmp(namesclean, 'CameraFrame'))
            Task_TTL_times = times(strcmp(namesclean, 'CameraFrame'));
        else
            Task_TTL_times = times(strcmp(namesclean, 'BeforeStimuliShown'));
        end
    else
        Task_TTL_times = times(strcmp(namesclean,'TTL'));
    end
    
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
    end
    % Then clip the start if it's bad
    badAO = find(diff(origAO) > 15, 1, 'last');
    if badAO < length(origAO)/2
        origAO = origAO(badAO+1:end);
    end
    newAO = origAO;
    
    % If there are still more received than sent task TTLs, try to scan
    % for them
    ldif = length(newAO) - length(Task_TTL_times);
    if ldif > 0
        ntry = ldif+1;
        trycorr = nan(ntry,1);
        for tryi = 1:ntry
            trycorr(tryi) = corr(diff(newAO((1:length(Task_TTL_times))+tryi-1))', diff(Task_TTL_times));
        end
        [~,maxi] = max(trycorr);
        newAO = newAO((1:length(Task_TTL_times))+maxi-1);
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
    times_conv = times*b(1) + b(2);
    
    % Store the task timestamps
    allevents = unique(namesclean);
    
    if any(strcmp(allevents, 'CameraFrame'))
        % Remove CameraFrame to do separately
        allevents(strcmp(allevents, 'CameraFrame')) = [];
        % Store CameraFrame info
        T.Frame_times = [nameidx(strcmp(namesclean,'CameraFrame')), times(strcmp(namesclean,'CameraFrame'))];
    end
    
    for k = 1:length(allevents)
        T.Time.(allevents{k}) = times_conv(strcmp(namesclean,allevents{k}));
    end
    
    % Get the target info
    T.Target_idx = nameidx(strcmp(namesclean, 'StartOfMovement'));
    
    angs_in_order = (0:7)*45;
    
    T.Target_dir = angs_in_order(T.Target_idx);
    T.Target_angle_zero = 'Right';
    T.Target_angle_orientation = 'Counter-Clockwise';
    
    % Load in kinematics
    % Get kinematic file name convention
    kinfiles = dir([Kindir '/*.csv']);
    kinsplit = strsplit(kinfiles(1).name,'session');
    datestr = kinsplit{1}(1:8);
    
    DLC = readcell([Kindir '/' kinsplit{1} 'session' num2str(Key.Vid(fii),'%0.3i') kinsplit{2}(4:end)]);

    if ndim == 2
        K.Labels = DLC(2,2:3:end);
        nlab = length(K.Labels);
        nframe = max(cell2mat(DLC(4:end,1)))+1;

        K.Xpos = nan(nframe, nlab);
        K.Ypos = nan(nframe, nlab);
        K.Likelihood = nan(nframe, nlab);
        for labi = 1:nlab
            DLCdat = DLC(4:end,[2 3 4]+((labi-1)*3));
            tempmask = cellfun(@ismissing,DLCdat);
            DLCdat(tempmask) = {nan};
            K.Xpos(:,labi) = cell2mat(DLCdat(:,1));
            K.Ypos(:,labi) = cell2mat(DLCdat(:,2));
            K.Likelihood(:,labi) = cell2mat(DLCdat(:,3));
        end
    elseif ndim == 3
        dlcheaders = DLC(1,:);
        for headi = 1:length(dlcheaders)
            tempsplit = strsplit(dlcheaders{headi},'_');
            if ~any(strcmp(tempsplit{1},{'M', 'center', 'fnum'}))
                stripheads{headi} = tempsplit{1};
            end
        end
        K.Labels = stripheads(1:6:end);
        nlab = length(K.Labels);
        nframe = max(cell2mat(DLC(2:end,end)))+1;
        K.Xpos = nan(nframe, nlab);
        K.Ypos = nan(nframe, nlab);
        K.Zpos = nan(nframe, nlab);
        K.Error = nan(nframe, nlab);
        K.Likelihood = nan(nframe, nlab);
        for labi = 1:nlab
            DLCdat = DLC(2:end,[1 2 3 4 6]+((labi-1)*6));
            tempmask = cellfun(@ismissing,DLCdat);
            DLCdat(tempmask) = {nan};
            K.Xpos(:,labi) = cell2mat(DLCdat(:,1));
            K.Ypos(:,labi) = cell2mat(DLCdat(:,2));
            K.Zpos(:,labi) = cell2mat(DLCdat(:,3));
            K.Error(:,labi) = cell2mat(DLCdat(:,4));
            K.Likelihood(:,labi) = cell2mat(DLCdat(:,5));
        end
    end

    % Rectify timestamps for kinematics...
    kintimes = cell2mat(readcell([Kintimedir '/session' num2str(Key.Vid(fii),'%0.3i') '/' datestr '_' SSi '_session' num2str(Key.Vid(fii),'%0.3i') '_topCam_timestamps.txt']));
    kintimes = kintimes(1:2:end)/(1e9);
    
    pauseAlign = 0;
    
    camAligni = find(diff(camTime)>1,1,'first');
    kinAligni = find(kintimes>1,1,'first');
    
    % Apply fixes to alignment indices
    if ~isempty(subfii)
        if ~isnan(fixtab{subfii,'pauseAlign'})
            pauseAlign = fixtab{subfii,'pauseAlign'};
            disp('Doing pause alignment');
        else
            if ~isnan(fixtab{subfii,'camAligni_direct'})
                camAligni = fixtab{subfii,'camAligni_direct'};
                disp('Did cam align index fix');
            end
            if ~isnan(fixtab{subfii,'camAligni_gap'})
                gaplocs = find(diff(camTime)>1);
                camAligni = gaplocs(fixtab{subfii,'camAligni_gap'}) + fixtab{subfii,'camAligni_gap_adj'};
                disp('Did cam align index fix');
            end
            if ~isnan(fixtab{subfii,'kinAligni'})
                kinAligni = fixtab{subfii,'kinAligni'};
                disp('Did kin align index fix');
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
    
    %         plot(K.Time,1:length(K.Time), '.r')
    if ~exist(savedir)
        mkdir(savedir)
    end
    
    save([savedir '/' datestr '_' SSi '_' num2str(fii,'%0.2i') '_Full_' num2str(ndim) 'D.mat'], 'N', 'T', 'K', '-v7.3');
end