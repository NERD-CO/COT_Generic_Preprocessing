% Run preprocessing steps once all data is generated and screened!

function Preprocess_All(datadir, regen, ndim, doplot, Settings)
    % Suppress the spline nan warning.
    warning('off','MATLAB:chckxy:IgnoreNaN');
    
    dddir = getenv('DBSDATADIR');

    regenall = any(strcmp(regen,'all'));
    
    % % Set and retrieve project settings
    % SetSS_COT_STN(datadir);
    % Settings = load([datadir 'SS.mat']);
    
    for subi = 1:length(Settings.SS)
        SSi = Settings.SS{subi};
    
        % Find the data tartares
        tardir = [dddir SSi '\Data_Tartare_Full_' num2str(ndim) 'D\'];
    
        if ~isfolder(tardir)
            mkdir(tardir);
        end
    
        tarfiles = dir([tardir '*' SSi '*Full_' num2str(ndim) 'D.mat']);
    
        % Make the data tartares if they don't yet exist
        if isempty(tarfiles) || any(strcmp(regen,'tar')) || regenall
            disp('');
            disp('Creating data tartares');
            create_data_tartare(dddir, SSi, tardir, [dddir 'COT_CamAlign_Fixes.xlsx'], ndim);
        else
            disp('Not regenerating data tartares');
        end
    
        tarfiles = dir([tardir '*' SSi '*Full.mat']);
    
        % The video annotation spreadsheets location
        annodir = [dddir SSi '\Annotated_Video\'];
    
        % Load the task - video - recording depth Key file for this subject
        keyfilename = [dddir SSi '\' SSi '_Key'];
        Key = tdfread(keyfilename);
    
        for fii = 1:length(tarfiles)
    
            savestr = [datadir 'Reach_Kin_Neu_' SSi '_' num2str(fii,'%0.2i') '.mat'];
    
            if ~exist(savestr, 'file') || any(strcmp(regen,'rkn'))
                fprintf('\n');
                disp(['Processing ' SSi ' file ' num2str(fii)]);
    
                % Load the data tartare
                Tar = load([tardir tarfiles(fii).name]);
                tarsplit = strsplit(tarfiles(fii).name,'_');
                tarnumstr = tarsplit{3};
                tarnum = str2num(tarnumstr);
    
                % Extract position and time from tartare file, touch it up
                [gapktime, rawpos] = extract_dlc_time_and_pos(Tar.K, Settings.kinname{subi}, Settings.likecutoff, ndim);
    
                % Fix skipped video frames
                [ktime, rawpos2] = fix_frameskips(gapktime, rawpos, Settings.vidframerate);
    
                % Fix blips and gaps in the kinematics
                [Pos] = fix_blips_and_gaps(ktime, rawpos2, Settings.vidframerate, 'endtime', Tar.T.Time.CenterIsShown(end));
    
                % Calculate higher order kinematics
                disp('');
                disp('Getting kinematics');
                [Time, Kin] = calculate_kinematics(ktime, Pos);
    
                % Get firing rates in ktime
                if isfield(Tar.N, 'SpkID') && ~isempty(Tar.N.SpkID)
                    disp('');
                    disp('Getting FRs');
                    [FR] = calculate_ktime_rates(Time.ktime, Tar.N, Settings.slidepad);
                else
                    FR = [];
                end
                %%
                % Attempt to detect reach starts and stops
                annofile = dir([annodir '\*_session' num2str(Key.Vid(tarnum),'%0.3i') '_*.xlsx']);
    
                % HIGHLY RECOMMENDED TO RUN THIS WITH DOPLOT FLAG AS 1, CLICK THROUGH THE PLOTS TO LOOK
                % FOR WEIRD THINGS HAPPENING, ADD X'S TO ANNOTATED VIDEO FILE
                % TO SCREEN WEIRD REACHES, OR EDIT COT_ReachFind_Fixes.xlsx AND
                % LATER YOU CAN RERUN IT WITH 0 ONCE YOU'VE SCREENED
                % EVERYTHING, if you have to regenerate (runs faster without
                % plotting)
                [Reach] = detect_reach_events(SSi, fii, Time.ktime, Kin.SmoothG50.Vel, Tar.T, Settings.vidframerate, [annodir annofile(1).name], [dddir 'COT_ReachFind_Fixes.xlsx'], doplot);
    
                % Save current stuff as a reach_kin_neu file
                disp('');
                N = Tar.N;

                % Also save task events from the tar file
                tosave = {'ArrowIsShown', 'CenterIsShown', 'StartOfMovement'};
                for tsi = 1:length(tosave)
                    Time.(tosave{tsi}) = Tar.T.Time.(tosave{tsi});
                end
                disp(['Saving ' savestr]);
                save(savestr, 'Reach', 'Time', 'Kin', 'FR', 'N');
                close all;
            end
        end
    end
end

%% Create the data tartare (mid processed data)
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
%
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
    
    % Where are the kinematic .csv files?
    % Kindir = ['C:/Users/RexBoxOne/Documents/DeepLabCut/dbsdlc/' SSi '/COTask-Rex-2021-07-07/videos'];
    if ndim == 2
        Kindir = dir([dddir '/' SSi '/' SSi '-*-*-*-*']);
        Kindir = [dddir '/' SSi '/' Kindir(1).name '/videos/'];
    elseif ndim == 3
        Kindir = [dddir '/' SSi '/3D_CSV/'];
    end
    
    % Where is the Key file?
    Keyfilename = [dddir '/' SSi '/' SSi '_Key'];
    
    % Read in the Key file
    Key = tdfread(Keyfilename);
    
    nfiles = length(Key.AODepth);
    
    for fii = 1:nfiles
        disp(['Doing file ' num2str(fii) ' (Depth ' num2str(Key.AODepth(fii)) ')']);
        
        % Load in the AO mat data
        AOfileload = dir([AOmatdir '/*D' num2str(Key.AODepth(fii), '%0.3f') 'F' num2str(Key.F(fii),'%0.4i') '*.mat']);
        M = load([AOmatdir '/' AOfileload.name]);
        
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
    
        [K.Time, T] = Align_TTLs_And_Frames(SSi, fii, Key, ndim);
        
        if ~exist(savedir)
            mkdir(savedir)
        end
        
        save([savedir '/' datestr '_' SSi '_' num2str(fii,'%0.2i') '_Full_' num2str(ndim) 'D.mat'], 'N', 'T', 'K', '-v7.3');
    end
end

%% Get things from DLC
function [gapktime, rawpos] = extract_dlc_time_and_pos(K, kinname, likecutoff, ndim)
    % Find the label in the DLC data
    labi = find(strcmp(K.Labels,kinname));
    if ndim == 2
        rawpos = [K.Xpos(:,labi), K.Ypos(:,labi)];
    elseif ndim == 3
        rawpos = [K.Xpos(:,labi), K.Ypos(:,labi), K.Zpos(:,labi)];
    end
    
    % snip the low likelihood
    rawpos(K.Likelihood(:,labi)<likecutoff,:) = NaN;
    disp(['Cut ' num2str(sum(K.Likelihood(:,labi)<likecutoff)) ' samples for low likelihood']);
    
    gapktime = K.Time;
    
    % Remove first two times to attempt to avoid start pauses
    gapktime = gapktime(3:end);
    rawpos = rawpos(3:end,:);
    
    % Try to detect a 1s end pause and stop before it.
    dgkt = diff(gapktime);
    lastsamp = find(dgkt > 0.98 & dgkt < 1.02);
    if ~isempty(lastsamp)
        if lastsamp > length(gapktime)*0.9
            gapktime = gapktime(1:lastsamp);
            rawpos = rawpos(1:lastsamp,:);
        end
    end
    disp('Detected end pause, clipping.');
end

%% Detect and fix skipped frames
% Function to perform fixing for skipped frames in the DLC video recording.
%
% Rex Tien 12/14/23

% Frame skip fixing:
%
% If the time vector looks fine, don't alter it.
%
% 1. Detect when the difference between frame times is a multiple of the
% framerate
%
% 2. Linear interpolation between frame times on either side of skips to
% fill in the time vector so that in the end, times occur at approximately
% the correct framerate.
%
% 3. Insert NaNs into the position matrix to represent missed tracking
% during the skipped frames

% Required inputs:
% gapktime: original time vector from DLC alignment (in seconds)
% rawpos: original positions matrix from DLC output (NxD matrix. As many
%       positions can be sent in as desired, each dimension as a column)
% framerate: the intended framerate of the video recording

% Optional settings:
% verbose: bool - Print things along the way? (default = true)

% Outputs:
% ktime: Time vector expanded to accomodate skipped frames
% pos: Position matrix with nans inserted where skipped frames existed

function [ktime, pos] = fix_frameskips(gapktime, rawpos, framerate, varargin)
    
    % Default values:
    V.verbose = true;
    
    % Get and check optional arguments
    nargin = length(varargin);
    if nargin > 0
        for argi = 1:2:nargin
            if ~isa(varargin{argi}, 'char')
                error('Malformed optional inputs - please enter as (''name'', value)');
            else
                V.(varargin{argi}) = varargin{argi+1};
            end
        end
    end
    
    boolchecks = {'verbose'};
    for booli = 1:length(boolchecks)
        if ~isa(V.(boolchecks{booli}), 'logical')
            error([boolchecks{booli} ' must be a boolean value (true or false)']);
        end
    end
    
    if size(gapktime,1) ~= size(rawpos,1)
        error('input time and raw position must have the same number of samples (rows)');
    end
    
    if any(isnan(gapktime))
        error('NaN value(s) detected in input time')
    end
    
    ndim = size(rawpos,2);
    if ndim<=0
        error('position matrix appears to be empty');
    end
    
    framestep = 1/framerate;
    
    % Do the frame skip fix if necessary
    dkt = diff(gapktime);
    % check for weird short gaps
    if any(unique(dkt) < 0.9*framestep)
        error('detected abnormally short time between frames - check time vector');
    end
    % Detect abnormally long gaps that would suggest missing video frames
    if any(unique(dkt) > 1.1*framestep) || any(unique(dkt) < 0.9*framestep)
        if V.verbose
            disp('Detected skipped frame(s), doing fix.');
        end
        % Convert time to samples by rounding
        fakeusetime = [0; cumsum(round(dkt/framestep))];
        fakestep = 1;
        dut = diff(fakeusetime);
    
        % Try to calculate the number of samples there should be
        nsampfix = round((fakeusetime(end)-fakeusetime(1))/fakestep)+1;
        fixusetime = nan(nsampfix,1);
        ktime = nan(nsampfix,1);
        addedpt = zeros(nsampfix,1);
        ndim = size(rawpos,2);
        pos = nan(nsampfix,ndim);
        fixcount = 1;
        for fixi = 1:length(dut)
            % if it's a good sample, store it and increment; if
            % it's bad store the current sample and the skipped
            % ones (not the next sample)
            nskipped = round(dut(fixi)/fakestep);
            % interpolate time during skipped frames
            addtime = interp1([fixi, fixi+nskipped], [fakeusetime(fixi), fakeusetime(fixi+1)], fixi:fixi+nskipped);
            fixusetime(fixcount:fixcount+nskipped-1) = addtime(1:end-1);
            addtimeorig = interp1([fixi, fixi+nskipped], [gapktime(fixi), gapktime(fixi+1)], fixi:fixi+nskipped);
            ktime(fixcount:fixcount+nskipped-1) = addtimeorig(1:end-1);
            % fill in nans in kinematics
            pos(fixcount:fixcount+nskipped-1,:) = [rawpos(fixi,:); nan(nskipped-1,ndim)];
            if (nskipped>1)
                addedpt(fixcount+1:fixcount+nskipped-1) = ones(nskipped-1,1);
            end
            fixcount = fixcount+nskipped;
        end
        fixusetime(end) = fakeusetime(end);
        ktime(end) = gapktime(end);
        pos(end,:) = rawpos(end,:);
        if V.verbose
            disp(['Fixed ' num2str(sum(addedpt)) ' missing frames']);
        end
    else
        % If no gaps, let it through
        if V.verbose
            disp('No skipped frames detected, returning');
        end
        ktime = gapktime;
        pos = rawpos;
    end
end

%% Function to process the raw position data by fixing blips and gaps and
% smoothing it. Data should be run through fix_frameskisp.m first to repair
% skipped frames. Function also smooths the position data (smoothed data
% should be used, as spline fill will be more effective on smoothed data!)
%
% Rex Tien 12/14/23

% Blip and Gap fixing:
%
% Definitions:
% Blips = short appearances of marker in otherwise bad tracking periods
% Gaps = short gaps due to missed frames or transiently low DLC confidence,
% or longer gaps due to DLC losing tracking completely
%
% Strategy:
% 1. Delete blips shorter than 'blipcut' UNLESS they are between two gaps
% that are shorter than 'gapfill'
% 1a. Special case: delete short blips that are between two short gaps if
% they are part of a "blip train" i.e. when several blips happen back to
% back.
%
% 2. Cut out any isolated tracking segments that occur more than 'afterend'
% after the last reach stop time. If afterend is 'NaN', don't perform this
% cutting.
%
% 3. Extend gaps longer than 'longthresh' by 'longcut' on either side.
%
% 4. Smooth the position data using Gaussian smoothing with kernel std
% deviations 'gkerns'
%
% 5. Fill gaps shorter than 'gapfill' in raw and smoothed positions using
% spline fill

% Required inputs:
% ktime: time vector with no frame skips (in seconds)
% rawpos: raw position vector after frame skip fixing (Nx2 or Nx3 matrix.
%       Each keypoint's position matrix must be sent separately to this function)
% framerate: the framerate of the video recording

% Optional settings:
% verbose: bool - Print things along the way? (default = true)
% doblips: bool - Do short blip deletion? (default = true)
% doshortgaps: bool - Delete short gaps? (default = true)
% dolonggaps: bool - Extend long gaps? (default = true)
% blipcut: Max length of short blips to cut (in seconds; default = 0.1)
% gapfill: Max length of short gaps to fill (in seconds; default = 0.1)
% endtime: The last behavioral timestamp. Position segments after this time
%       will be cut (no default - this value must be supplied if you want
%       to do end cutting)
% afterend: The the after 'endtime' to start cutting position segments (in
%       seconds; default = 3)
% longthresh: Min length of long gaps to extend (in seconds; default = 5)
% longcut: Amount to extend long gaps by on either side (in seconds; default = 0.5)
% gkerns: Std deviations of the Gaussian kernel for smoothing, can be
%       arbitrary length vector, will run smoothing multiple times with
%       each kernel setting.

% Outputs:
% Pos, struct with fields:
%   Raw (blip and gap fixed raw positions (i.e. not smoothed))
%   SmoothGXX (smooth blip and gap fixed positions, with a field for each
%       input gkern, where XX is num2str of gkern, converted to milliseconds (gkern of 0.05 becomes G50))

function [Pos] = fix_blips_and_gaps(ktime, pos, framerate, varargin)
    
    % Default values:
    V.verbose = true;
    V.doblips = true;
    V.doshortgaps = true;
    V.dolonggaps = true;
    V.blipcut = 0.1;
    V.gapfill = 0.1;
    V.afterend = 3;
    V.longthresh = 5;
    V.longcut = 0.5;
    V.gkerns = [0.015, 0.05];
    
    % Get and check optional arguments
    nargin = length(varargin);
    if nargin > 0
        for argi = 1:2:nargin
            if ~isa(varargin{argi}, 'char')
                error('Malformed optional inputs - please enter as (''name'', value)');
            else
                V.(varargin{argi}) = varargin{argi+1};
            end
        end
    end
    
    boolchecks = {'verbose', 'doblips', 'doshortgaps', 'dolonggaps'};
    for booli = 1:length(boolchecks)
        if ~isa(V.(boolchecks{booli}), 'logical')
            error([boolchecks{booli} ' must be a boolean value (true or false)']);
        end
    end
    
    numchecks = {'blipcut', 'gapfill', 'afterend', 'longthresh', 'longcut', 'gkerns'};
    for numi = 1:length(numchecks)
        if ~isa(V.(numchecks{numi}), 'numeric')
            error([numchecks{numi} ' must be a number']);
        elseif ~(V.(numchecks{numi}) > 0)
            error([numchecks{numi} ' must be greater than 0']);
        end
    end
    
    if isfield(V,'endtime')
        if ~isa(V.endtime, 'numeric')
            error('endtime must be a number');
            if ~(V.endtime > 0)
                error('endtime must be greater than 0');
            end
        end
    end
    
    if size(ktime,1) ~= size(pos,1)
        error('input time and position must have the same number of samples (rows)');
    end
    
    if any(isnan(ktime))
        error('NaN value(s) detected in input time');
    end
    
    ndim = size(pos,2);
    if (ndim<=0) || (ndim>3)
        error('position matrix must be Nx1, Nx2 or Nx3');
    end

    if ~all(any(isnan(pos),2) == all(isnan(pos),2))
        error(['Detected different gaps in the different dimensions of input' ...
            'position matrix. This may indicate that not all positions are from the same DLC keypoint - please check!']);
    end
    
    framestep = 1/framerate;
    if any(unique(diff(ktime)) > 1.1*framestep)
        error('dropped frames detected in the time vector, please run fix_frameskips.m first');
    end
    
    blipcut_frames = V.blipcut*framerate;
    gapfill_frames = V.gapfill*framerate;
    longthresh_frames = V.longthresh*framerate;
    longcut_frames = V.longcut*framerate;
    
    if V.doblips
        % Find blips
        [blipstarts, blipstops, bliplens] = findblipgaps(pos,1);
        cuts = find(bliplens < blipcut_frames);
        ncut = length(cuts);
        nblips1 = sum(bliplens < framerate);
        
        % Find gaps
        [gapstarts, gapstops, gaplens] = findblipgaps(pos, 0);
        fills = find(gaplens < gapfill_frames);
        smallgapstarts = gapstarts(fills);
        smallgapstops = gapstops(fills);
        ngaps1 = length(gaplens);
        
        % Cut out the blips. if a blip is between two short gaps, don't delete
        % it, unless it's in a series of blipgaps (a bliptrain)!
        nspared = 0;
        nactuallycut = 0;
        for blipi = 1:ncut
            docut = 1;
            beforestart = blipstarts(cuts(blipi))-1;
            afterstop = blipstops(cuts(blipi))+1;
            % Detect if it's between two small gaps
            if (any(smallgapstops==beforestart) && any(smallgapstarts==afterstop))
                prevgapstart = smallgapstarts(smallgapstops==beforestart);
                nextgapstop = smallgapstops(smallgapstarts==afterstop);
                % Detect if there's another blip down the line
                if (blipi ~= ncut) && (nextgapstop == blipstarts(cuts(blipi+1))-1)
                    docut = 1;
                % Detect if there's another blip up the line
                elseif (blipi > 1) && (prevgapstart == blipstops(cuts(blipi-1))+1)
                    docut = 1;
                else
                    docut = 0;
                end
            end
        
            if docut
                pos(blipstarts(cuts(blipi)):blipstops(cuts(blipi)),:) = NaN;
                nactuallycut = nactuallycut+1;
            else
                nspared = nspared+1;
            end
        end
        
        if V.verbose
            disp(['Dectected ' num2str(sum(bliplens < 1000)) ' blips < 1 second long, ' num2str(ncut) ' blips < ' num2str(V.blipcut) ' seconds.']);
            if ncut > 0
                disp(['Cut ' num2str(nactuallycut) ' blips shorter than ' num2str(V.blipcut) ' seconds.']);
                disp(['Spared ' num2str(nspared) ' blips that were between two short gaps.']);
            else
                disp('No blips were cut.');
            end
        end
    else
        if V.verbose
            disp('Not doing blip deletion.');
        end
    end
    
    if isfield(V,'endtime')
        % Cut out any blips that happen at least 3 seconds after the last center
        % show
        endcutoff = find(ktime > V.endtime + V.afterend,1,'first');
        if ~isempty(endcutoff)
            endblips = find(blipstarts > endcutoff);
            if ~isempty(endblips)
                pos(blipstarts(endblips(1)):end,:) = NaN;
                if V.verbose
                    disp(['Cut ' num2str(length(endblips)) ' segments at the end of the recording.']);
                end
            end
        end
    else
        if V.verbose
            disp('Not doing end cuts.');
        end
    end
    
    if V.dolonggaps
        % Now redo the gap finding for long gap extension
        [gapstarts, gapstops, gaplens] = findblipgaps(pos, 0);
    
        % Trim ends that are adjacent to very long gaps
        longs = find(gaplens > longthresh_frames); % Find gaps longer than 1 second
        nlong = length(longs);
        if V.verbose
            disp(['Detected ' num2str(nlong) ' long gaps longer than ' num2str(V.longthresh) ' seconds']);
        end
        for longi = 1:nlong
            if gapstarts(longs(longi)) > longcut_frames
                pos(gapstarts(longs(longi))-longcut_frames:gapstarts(longs(longi))-1,:) = NaN;
            else
                pos(1:gapstarts(longs(longi))-1,:) = NaN;
            end
            if gapstops(longs(longi)) < length(ktime)-longcut_frames
                pos(gapstops(longs(longi))+1:gapstops(longs(longi))+longcut_frames,:) = NaN;
            else
                pos(gapstops(longs(longi))+1:end,:) = NaN;
            end
        end
        if V.verbose
            disp(['Extended ' num2str(nlong) ' long gaps by ' num2str(V.longcut) ' seconds on each side']);
        end
    else
        if V.verbose
            disp('Not doing long gap extension.');
        end
    end
    
    % Now smooth. Retain the nan locations. Do it for every entry in gkerns
    nkerns = length(V.gkerns);
    kns = cell(nkerns,1);
    for kerni = 1:nkerns
        kns{kerni} = ['SmoothG' num2str(V.gkerns(kerni)*1000)];
        Pos.(kns{kerni}) = GaussSmooth_Arbitrary(ktime, pos, ktime, V.gkerns(kerni), 5);
        Pos.(kns{kerni})(isnan(pos(:,1)),:) = NaN;
    end
    
    if V.doshortgaps
        % Now redo the gap finding for short gap filling
        [gapstarts, gapstops, gaplens] = findblipgaps(pos, 0);
        fills = find(gaplens < gapfill_frames);
        nfill = length(fills);
        
        % Then, fill gaps with spline if they are < gapfill and isolated
        fillis = [];
        for gapi = 1:nfill
            fillis = [fillis, gapstarts(fills(gapi)):gapstops(fills(gapi))];
        end
        
        if any(isnan(pos(:,1)) ~= isnan(Pos.(kns{1})(:,1)))
            error('different gaps in pos and smoothpos for some reason');
        end

        nni = ~isnan(pos(:,1));
        nnt = ktime(nni);
        nnpos = pos(nni,:);
        for dimi = 1:ndim
            pos(fillis,dimi) = spline(nnt, nnpos(:,dimi), ktime(fillis));
        end
        for kerni = 1:nkerns
            nnsmoothpos = Pos.(kns{kerni})(nni,:);
            for dimi = 1:ndim
                Pos.(kns{kerni})(fillis,dimi) = spline(nnt, nnsmoothpos(:,dimi), ktime(fillis));
            end
        end
        
        % detect if filling somehow failed
        if any(any(isnan(pos(fillis,:)))) || any(any(isnan(Pos.(kns{1})(fillis,:))))
            error('Somehow gapfilling failed!')
        end
        if V.verbose
            disp(['Filled ' num2str(nfill) ' short gaps shorter than ' num2str(V.gapfill) ' seconds.']);
        end
    else
        if V.verbose
            disp('Not doing short gap filling.');
        end
    end
        
    % Detect blips and gaps agian one more time for fun
    [~, ~, remainingbliplens] = findblipgaps(pos, 1);
    [~, ~, remaininggaplens] = findblipgaps(pos, 0);
    % If there are somehow sitll gaps left
    if any(remaininggaplens < gapfill_frames)
        error('Somehow still small gaps after gap filling!')
    end

    if V.verbose
        disp(['All done! ' num2str(length(remaininggaplens)) ' gaps remain.']);
        disp([num2str(sum(remainingbliplens<framerate)) ' blips < 1 second remain.']);
        disp([num2str(length(remaininggaplens)) ' gaps remain.']);
    end

    % Make sure to store rawpos too!!!!!!
    Pos.Raw = pos;
end

function [starts, stops, lengths] = findblipgaps(pos, blipgapflag)
    % If blipgapflag is 1, find blips. If 0 find gaps
    if blipgapflag == 1
        locs = ~isnan(pos(:,1));
    elseif blipgapflag == 0
        locs = isnan(pos(:,1));
    else
        error('Bad blipgapflag');
    end

    dl = diff(locs);
    starts = find(dl==1)+1;
    if locs(1) == 1
        starts = [1; starts];
    end
    stops = find(dl==-1);
    if locs(end) == 1
        stops = [stops; length(locs)];
    end

    lengths = stops-starts + 1;
end

%%
% Function to get higher order kinematics (velocity, acceleration and
% associated terms) from smooth position input.
%
% Rex Tien 12/14/23

% Kinematic processing
%
% Differentiation is done using [diff] function
%
% Smoothing is done using GaussSmooth_Arbitrary and is done zero times. The 
% smoothed position is used as input to generate 'smooth' higher order kinematics, 
% and raw is used to generate raw. Higher order kinematics are resampled
% back to position time using [spline].

% Required inputs:
% ktime: Aligned timestamps from DLC after frame skip fixing
% Pos: struct with fields:
%       Raw: frame-fixed positions matrix, not smoothed
%       Smooth.GXX: frame-fixed positions matrix, smoothed with Gaussian
%           kernel of width XX milliseconds, can have any number of GXX fields

% Outputs:
% Time: Struct with fields
%   ktime: same as input ktime
%
% Kin: Struct with kinematics resampled to position times (ktime)
%   Raw: Non-smoothed - recommended not to use as higher order terms will be very jumpy
%       pos
%       vel
%       speed
%       acc
%       accmag: magnitude of acceleration
%       accsign: magnitude of acceleration multiplied by the dot product between velocity and acceleration (positive when increasing the current acceleration, negative when decreasing it)
%   SmoothGXX: Kinematics calculated from smooth position (Smoothed with Gauss kernel of width XX milliseconds) with no re-smoothing
%       (same subfields as Kin.Raw)

function [Time, Kin] = calculate_kinematics(ktime, Pos)

    % Do some input checking
    
    if size(ktime,1) ~= size(Pos.Raw,1)
        error('input time and position must have the same number of samples (rows)');
    end
    
    if any(isnan(ktime))
        error('NaN value(s) detected in input time')
    end
    
    ndim = size(Pos.Raw,2);
    if ndim<=0
        error('position matrix appears to be empty');
    end
    
    smoothtypes = fieldnames(Pos);
    
    % Setting up
    ndim = size(Pos.Raw,2);
    Time.ktime = ktime;
    vtime = mean([ktime(1:end-1), ktime(2:end)], 2);
    atime = mean([vtime(1:end-1), vtime(2:end)], 2);

    % Get velocity, acc, resample back to ktime, calculate other kinematics
    for sti = 1:length(smoothtypes)
        stype = smoothtypes{sti};

        Kin.(stype).Pos = Pos.(stype);

        tempvel = diff(Kin.(stype).Pos,1)./diff(ktime);
        tempacc = diff(tempvel,1)./diff(vtime);

        Kin.(stype).Vel = nan(length(ktime),ndim);
        for dimi = 1:ndim
            Kin.(stype).Vel(:,dimi) = spline(vtime, tempvel(:,dimi), ktime);
        end
        Kin.(stype).Vel(isnan(Kin.(stype).Pos(:,1)),:) = NaN;
        Kin.(stype).Speed = vecnorm(Kin.(stype).Vel')';
        Kin.(stype).Dirs = Kin.(stype).Vel./(repmat(Kin.(stype).Speed,[1,ndim]));

        Kin.(stype).Acc = nan(length(ktime),ndim);
        for dimi = 1:ndim
            Kin.(stype).Acc(:,dimi) = spline(atime, tempacc(:,dimi), ktime);
        end
        Kin.(stype).Acc(isnan(Kin.(stype).Pos(:,1)),:) = NaN;
        Kin.(stype).AMag = vecnorm(Kin.(stype).Acc')';

        avdot = sum((Kin.(stype).Vel./Kin.(stype).Speed).*(Kin.(stype).Acc./Kin.(stype).AMag),2);
        Kin.(stype).ASign = Kin.(stype).AMag.*avdot;
    end
end

%%
function [FR] = calculate_ktime_rates(ktime, N, slidepad)
    V.gkerns = [0.015, 0.05];
    
    nneu = length(N.SpkTimes);
    nkp = length(ktime);
    
    % Do ms resolution FIR rates, add a 5 second buffer on the front and back
    % to avoid edge effects when smoothing
    msedge = ((ktime(1)-5):0.001:(ktime(end)+5))';
    mslooktime = msedge(1:end-1)+0.0005;
    
    msrates = nan(length(mslooktime), nneu);
    for j = 1:nneu
        msrates(:,j) = Spikes2FIR_Arbitrary( N.SpkTimes{j}, msedge)';
    end
    
    for kerni = 1:length(V.gkerns)
        thisname = ['SmoothG' num2str(V.gkerns(kerni)*1000,'%0.2i')];
        FR.(thisname).krates = nan(nkp,nneu);
        for neui = 1:nneu
            FR.(thisname).krates(:,neui) = GaussSmooth_Arbitrary( mslooktime, msrates(:,neui), ktime, V.gkerns(kerni), 5);
        end
    
        padrates = [nan(slidepad,nneu); FR.(thisname).krates; nan(slidepad,nneu)];
    
        if any(any(isnan(FR.(thisname).krates)))
            OHNOO
        end
    
        % And do sliding renormalized and mean subtracted rates
        FR.(thisname).slidenormrates = nan(size(FR.(thisname).krates));
        FR.(thisname).slidemeanrates = nan(size(FR.(thisname).krates));
        for slidei = 1:nkp
            FR.(thisname).slidemeanrates(slidei,:) = (FR.(thisname).krates(slidei,:) - nanmean(padrates(slidei:(slidei+slidepad*2),:),1));
            FR.(thisname).slidenormrates(slidei,:) = FR.(thisname).slidemeanrates(slidei,:)./nanstd(padrates(slidei:(slidei+slidepad*2),:),0,1);
        end
    end
end

%%
% Attempt to detect reach events (reach starts, stops, peak times)
% Requires:
%   ktime - frame fixed kinematic timestamps
%   vel - smoothed, frame fixed velocity
%   T - tartare time struct with behavioral event times

function [Reach] = detect_reach_events(SSi, fii, ktime, vel, T, framerate, annofile, fixfile, doplot)
    
    ndim = size(vel,2);

    plotpeakfindin = doplot;
    plotpeakfindout = doplot;
    plotstartstopfind = doplot;

    nout = length(T.Time.StartOfMovement);
    nback = length(T.Time.CenterIsShown);

    % Load video notes and find x'd out reaches
    vidnotes = readcell(annofile);
    if size(vidnotes,2) < 4 || ~strcmp(vidnotes{1,4}, 'g')
        error('Annotated Video has not been reviewed (no ''g'' in cell 1,4)');
    end
    xout = strcmp(vidnotes(strcmp(vidnotes(:,2),'Out'),3),'x');
    xback = strcmp(vidnotes(strcmp(vidnotes(:,2),'Back'),3),'x');
    lout = strcmp(vidnotes(strcmp(vidnotes(:,2),'Out'),3),'l');
    lback = strcmp(vidnotes(strcmp(vidnotes(:,2),'Back'),3),'l');

    % Load fixes file and extract relevant rows
    fixnotes = readcell(fixfile);
    fixheaders = fixnotes(1,:);
    fixtab = cell2table(fixnotes(2:end,:));
    fixtab.Properties.VariableNames = fixheaders;
    subfii = strcmp(fixtab{:,'Subject'},SSi) & fixtab{:,'Filenum'}==fii;

    % Low-pass filter the velocity and recalculate speed before doing
    % event hunting to try to avoid the effects of tremor.
    filtcutoff = 3;
    [filtB, filtA] = butter(4, filtcutoff/(framerate/2), 'low');
    % De-nan the vel
    interpvel = vel;
    fills = isnan(vel(:,1));
    for dimi = 1:size(vel,2) % Interpolate the nans so smoother can work
        interpvel(fills,dimi) = interp1(ktime(~fills), vel(~fills,dimi), ktime(fills));
    end
    interpvel(isnan(interpvel)) = 0;
    filtvel = filtfilt(filtB, filtA, interpvel);
    filtvel(isnan(vel(:,1)),:) = nan; % Put back the  nans
    filtspeed = vecnorm(filtvel')';
    filtdirs = filtvel./repmat(filtspeed,[1,ndim]);

    reachpeaks_rough = [];
    reachnum = [];
    reachdirs = [];
    outreach = [];
    reactiontime = [];

    % Find peaks in -1 to 1.5s around center show time
    for i = 1:nback
        if ~xback(i)
            cback = 0.5;
            cfwd = 2.5;
            if lback(i)
                cback = 0;
            end

            fixcback = fixtab{subfii & fixtab{:,'reachnum_peaks'}==i, 'cback'};
            if ~isempty(fixcback) && ~isnan(fixcback)
                cback = fixcback;
            end
            fixcfwd = fixtab{subfii & fixtab{:,'reachnum_peaks'}==i, 'cfwd'};
            if ~isempty(fixcfwd) && ~isnan(fixcfwd)
                cfwd = fixcfwd;
            end

            ctime = T.Time.CenterIsShown(i);
            tmask = (ktime > ctime-cback) & (ktime < ctime+cfwd);
            [themax,maxi] = max(filtspeed(tmask));
            thisktime = ktime(tmask);
            if ~isempty(thisktime(maxi)) && any(~isnan(filtspeed(tmask)))
                reachpeaks_rough = [reachpeaks_rough; thisktime(maxi)];
                reachnum = [reachnum; i];
                reachdirs = [reachdirs; T.Target_dir(i)-180];
                outreach = [outreach; 0];

                if plotpeakfindin
                    figure(1);
                    clf;
                    hold on;
                    plot(thisktime, filtspeed(tmask));
                    plot(ctime, 1, '*b');
                    plot(thisktime(maxi), themax, '*r');
                    title(['Reach #' num2str(reachnum(end)) ' Reach Back']);
                    while ~waitforbuttonpress
                    end
                end
            end
        end
    end

    % Find peaks in -0.5 to 2.5s around target show time
    for i = 1:nout
        if ~xout(i)
            tback = 0.5;
            tfwd = 2.5;
            if lout(i)
                tback = 0;
            end

            fixtback = fixtab{subfii & fixtab{:,'reachnum_peaks'}==i, 'tback'};
            if ~isempty(fixtback) && ~isnan(fixtback)
                tback = fixtback;
            end
            fixtfwd = fixtab{subfii & fixtab{:,'reachnum_peaks'}==i, 'tfwd'};
            if ~isempty(fixtfwd) && ~isnan(fixtfwd)
                tfwd = fixtfwd;
            end

            ttime = T.Time.StartOfMovement(i);
            tmask = (ktime > ttime-tback) & (ktime < ttime+tfwd);
            [themax,maxi] = max(filtspeed(tmask));
            thisktime = ktime(tmask);
            if ~isempty(thisktime(maxi))
                reachpeaks_rough = [reachpeaks_rough; thisktime(maxi)];
                reachnum = [reachnum; i];
                reachdirs = [reachdirs; T.Target_dir(i)];
                outreach = [outreach; 1];

                if plotpeakfindout % && any(~isnan(filtspeed(tmask)))
                    figure(1);
                    clf;
                    hold on;
                    plot(thisktime, filtspeed(tmask));
                    plot(ttime, 1, '*b');
                    plot(thisktime(maxi), themax, '*r');
                    title(['Reach #' num2str(reachnum(end)) ' Reach Out']);
                    while ~waitforbuttonpress
                    end
                end
            end
        end
    end

    if any(isnan(reachpeaks_rough))
        error('NAN REACHPEAKS!');
    end

    [reachpeaks_rough, sorti] = sort(reachpeaks_rough);
    reachdirs = reachdirs(sorti);
    reachnum = reachnum(sorti);
    outreach = outreach(sorti);

    nreach = length(reachpeaks_rough);

    % Now that we have roughly found the reach peaks and put them in
    % order, use finer search to find reachstarts, reachstops and
    % reachpeak times
    reachstarts = nan(nreach,1);
    reachstops = nan(nreach,1);
    reachpeaks = nan(nreach,1);
    reactiontime = nan(nreach,1);

    % Look period is 2 second in front of and 2 s behind peaks
    lookback = 2;
    lookfwd = 2;
    reachdirwin = 0.5;
    reachdirsamps = reachdirwin*framerate/2;
    for k = 1:nreach
        peaki = find(ktime == reachpeaks_rough(k));
        backi = find((ktime >= reachpeaks_rough(k)-lookback) & (ktime <= reachpeaks_rough(k)));
        fwdi = find((ktime <= reachpeaks_rough(k)+lookfwd) & (ktime > reachpeaks_rough(k)));
        looki = [backi; fwdi];
        didproj = 0;

        mslooktime = ktime(looki(1)):0.001:ktime(looki(end));
        msbacktime = mslooktime(mslooktime <= reachpeaks_rough(k));
        msfwdtime = mslooktime(mslooktime > reachpeaks_rough(k));
        
        isallow = false;
        dontallow = false;

        nanplot = false;
        startstopplot = false;

        % if statement to exclude ones too close to end of recording
        if (peaki+reachdirsamps < length(filtvel)) && (peaki-reachdirsamps > 0)

            if any(fixtab{subfii,'reachnum_startstop'} == k)
                isallow = fixtab{subfii & fixtab{:,'reachnum_startstop'}==k, 'allownan'} == 1;
                dontallow = fixtab{subfii & fixtab{:,'reachnum_startstop'}==k, 'dontallow'} == 1;
            end

            % if statement to exclude ones where there is nan (Special
            % fixes for ones which actually seem fine)
            if ~any(isnan(filtvel(looki,1))) || isallow
                    
                didproj = 1;
                % 3D edit
                mv = mean(filtvel(peaki-reachdirsamps:peaki+reachdirsamps,:),1);
                
                reachdir = mv/vecnorm(mv);
                
                backv = filtspeed(backi).*(filtdirs(backi,:)*reachdir');
                fwdv = filtspeed(fwdi).*(filtdirs(fwdi,:)*reachdir');
                lookv = filtspeed(looki).*(filtdirs(looki,:)*reachdir');

                mslookv = interp1(ktime(looki), lookv, mslooktime);
                msbackv = mslookv(mslooktime <= reachpeaks_rough(k));
                msfwdv = mslookv(mslooktime > reachpeaks_rough(k));

                % Just pick out the zero crossings or a negative peak
                % that is <0.05 max, but now do it with the ms resampled
                % velocity!!!
                [maxrdspeed, maxrdi] = max(mslookv);

                reachpeaks(k) = mslooktime(maxrdi);

                if abs(reachpeaks(k) - reachpeaks_rough(k)) > 0.25
                    warning('Large offset between peak speed and peak speed in reach direction!')
                end

                if sign(backv(end))==-1
                    error('Negative reach speed near peak??');
                end

                version = ver('MATLAB');
                if contains(version.Release, '2022')
                    bxi = findpeaks(-abs(msbackv));
                    loci = find(msbackv(bxi.loc) < 0.05*maxrdspeed,1,'last');
                    if ~isempty(loci)
                        reachstarts(k) = msbacktime(bxi.loc(loci));
                    end

                    fxi = findpeaks(-abs(msfwdv));
                    loci = find(msfwdv(fxi.loc)<0.05*maxrdspeed,1,'first');
                    if ~isempty(loci)
                        reachstops(k) = msfwdtime(fxi.loc(loci));
                    end
                else
                    [~,bxi] = findpeaks(-abs(msbackv));
                    loci = find(msbackv(bxi) < 0.05*maxrdspeed,1,'last');
                    if ~isempty(loci)
                        reachstarts(k) = msbacktime(bxi(loci));
                    end

                    [~,fxi] = findpeaks(-abs(msfwdv));
                    loci = find(msfwdv(fxi)<0.05*maxrdspeed,1,'first');
                    if ~isempty(loci)
                        reachstops(k) = msfwdtime(fxi(loci));
                    end
                end

                if outreach(k)
                    reactiontime(k) = reachstarts(k) - T.Time.StartOfMovement(reachnum(k));
                else
                    reactiontime(k) = reachstarts(k) - T.Time.CenterIsShown(reachnum(k));
                end

                startstopplot = any(isnan([reachstarts(k) reachstops(k)]));
            else
                nanplot = true;
            end
            if (nanplot && ~dontallow) || startstopplot
                alli = ([backi; fwdi]);
                figure(1); clf; hold on;

                plot(ktime(alli), filtspeed(alli)); hold on;
                plot(ktime(ktime==reachstarts(k)), filtspeed(ktime==reachstarts(k)), '*r')
                plot(ktime(ktime==reachpeaks_rough(k)), filtspeed(ktime==reachpeaks_rough(k)), '*g')
                plot(ktime(ktime==reachstops(k)), filtspeed(ktime==reachstops(k)), '*b')

                plot([ktime(alli(1)) ktime(alli(1))], [0 100], '-k');
                plot([ktime(alli(end)) ktime(alli(end))], [0 100], '-k');

                plot(get(gca, 'xlim'), [0, 0], '-k');

                warnstr = [];
                if nanplot
                    warnstr = 'NAN DETECTED';
                elseif startstopplot
                    warnstr = 'COULDN''T FIND START OR STOP';
                end
                if outreach(k)
                    outstr = 'Out';
                else
                    outstr = 'In';
                end
                title([num2str(k) ' out of ' num2str(nreach) ' (Trial ' num2str(reachnum(k)) ' ' outstr ') ' warnstr]);

                while ~waitforbuttonpress
                end
            end
        end

        % % Step through all reaches in the 3D projected reach vel
        if plotstartstopfind && didproj
            alli = ([backi; fwdi]);
            allv = ([backv; fwdv]);
            allt = ktime(alli);

            figure(1); clf; hold on;
            if outreach(k)
                outstr = 'Out';
            else
                outstr = 'In';
            end
            if find(mslooktime == reachstops(k)) == length(mslooktime) || find(mslooktime == reachstarts(k)) == 1
                limstr = 'START OR STOP IS AT THE LIMIT!!';
            else
                limstr = [];
            end
            title([num2str(k) ' out of ' num2str(nreach) ' (Trial ' num2str(reachnum(k)) ' ' outstr ') ' limstr]);

            plot(allt, allv); hold on;

            plot(reachstarts(k), mslookv(mslooktime==reachstarts(k)), '*r')
            plot(reachpeaks(k), mslookv(mslooktime==reachpeaks(k)), '*g')
            plot(reachstops(k), mslookv(mslooktime==reachstops(k)), '*b')

            plot(mslooktime, mslookv);

            plot(get(gca, 'xlim'), [0, 0], '-k');


            if ~isempty(limstr)
                STOPPPPPPPPPP
            end
            while ~waitforbuttonpress
            end
        end
    end

    % Stop for investigation if there are nan reachstarts or stops
    
    if (any(isnan(reachstarts)) || any(isnan(reachstops(1:end-1)))) && ...
            ~any(fixtab{subfii, 'dontallow'}==1)
        badreachstarts = find(isnan(reachstarts))
        badreachstops = find(isnan(reachstops))
        error('Bad reachfind!')
    end

    % Clip the last reachstop if it's at the end of the recording :(
    if (reachstops(end) == ktime(end)) %|| isnan(reachstops(end))
        reachstops(end) = nan;
        disp('CLIPPED LAST REACHSTOP!');
    end

    % Clip reaches with NaN reach start or stop
    % badreach = find(isnan(reachstarts) | isnan(reachstops));
    badreach = find(isnan(reachstarts) | isnan(reachstops)); %| (reachstops-reachstarts < outlier_limits(subi,1)) | (reachstops-reachstarts > outlier_limits(subi,2))); % CHANGED THIS LINE FOR OUTLIER EXCLUSION
    reachdirs(badreach) = [];
    reachnum(badreach) = [];
    outreach(badreach) = [];
    reachpeaks_rough(badreach) = [];
    reachpeaks(badreach) = [];
    reachstarts(badreach) = [];
    reachstops(badreach) = [];
    reactiontime(badreach) = [];

    if any(groupcounts(reachstarts) > 1)
        error('Repeat reachstarts detected!');
    end

    store = {'reachdirs', 'reachnum', 'outreach', 'reachpeaks_rough', 'reachpeaks', 'reachstarts', 'reachstops', 'reactiontime'};
    for stori = 1:length(store)
        eval(['Reach.(store{stori}) = ' store{stori} ';']);
    end
end