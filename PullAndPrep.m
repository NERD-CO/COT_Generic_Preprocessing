% Initialize the subject data folder and do the first few steps of
% preprocessing!

% sstrs is a cell array of strings holding the subject codes to do (i.e.
% {'VM02', 'VM05', 'PD03', 'GP01'}

function PullAndPrep(sstrs, doplx, dodlc, ndim)
    for i = 1:length(sstrs)
        sstr = sstrs{i};

        PullFromSynology(sstr);
        AutoGenKey(sstr);
        if doplx
            disp('Converting .mat AO files to .plx for sorting');
            Batch_AO2plx(sstr);
        end
        if dodlc
            GenRunDLCPy(sstr, ndim);
        end
    end
end

% Pull files from synology, make nice folder, prep the alphaomega data

function PullFromSynology(subjstr)
    disp(['Pulling files for ' subjstr]);

    syndir = getenv('SYNOLOGYDIR');
    dddir = getenv('DBSDATADIR');
    aoconverterdir = getenv('AOCONVERTERDIR');

    abort = 0;

    % Check if folder exists
    if exist([dddir '\' subjstr], 'dir')
        disp('Folder already exists, not gonna mess with it');
    else
        mkdir([dddir '\' subjstr]);
        % Check if all the files are there and ready to go
        viddir = [syndir '\Video\' subjstr '_Compressed_Video'];
        if ~exist(viddir, 'dir')
            disp('Video files not found, aborting');
            abort = 1;
        end

        aodir = dir([syndir '\AO_Case_Data\' subjstr '*']);
        if isempty(aodir)
            disp('AO data not found, aborting');
            abort = 1;
        end

        taskdir = [syndir '\Task\' subjstr '_Task'];
        if ~exist(taskdir, 'dir')
            disp('Task files not found, aborting');
            abort = 1;
        end

        notefile = [syndir '\Notes\' subjstr '_Notes'];
        if ~exist(notefile, 'file')
            notefile = [syndir '\Notes\' subjstr '_Notes.ods'];
            if ~exist(notefile, 'file')
                disp('Notes file not found, aborting');
                abort = 1;
            end
%             
        end

        if ~abort
            % Pull the files
            copyfile(viddir, [dddir '\' subjstr '\' subjstr '_Compressed_Video']);
            copyfile([syndir '\AO_Case_Data\' aodir(1).name], [dddir '\' subjstr '\AO_Case_Data']);
            mkdir([dddir '\' subjstr '\AO_Case_Data\matfiles']);
            copyfile(taskdir, [dddir '\' subjstr '\' subjstr '_Task']);
            if exist(notefile, 'file')
                copyfile(notefile, [dddir '\' subjstr '\' subjstr '_Notes']);
            end
            disp('Done!');
        end
    end

    disp('Converting the AO .mpx files to .mat');
%     Convert the AO files
    mpxlist = dir([dddir '\' subjstr '\AO_Case_Data\*.mpx']);
    matdest = [dddir '\' subjstr '\AO_Case_Data\matfiles'];
    currentfold = pwd;
    cd(aoconverterdir);
    for i = 1:length(mpxlist)
        fprintf([num2str(i,'%0.3i') '/' num2str(length(mpxlist),'%0.3i')]);
        mpxsc = [dddir '\' subjstr '\AO_Case_Data\' mpxlist(i).name];
        command=['Mapfile.exe =' mpxsc '=Matlab=' matdest '='];
        [status,cmdout] = dos(command);
        fprintf('\b\b\b\b\b\b\b');
    end
    cd(currentfold);
    fprintf('\n');
    
    % Check for broken-up AO session files, fix them if found - added
    % 11/23/22 by Rex to accomodate Starfix cases which were split up
    matlist = dir([matdest '\*.mat']);
    nmat = length(matlist);
    mlist = cell(nmat,1);
    for fii = 1:nmat
        mlist{fii} = matlist(fii).name;
    end
    depthlist = extractBetween(mlist,'D','F');
    udep = unique(depthlist);
    ndep = length(udep);
    for depi = 1:ndep
        depthmatch = find(strcmp(depthlist, udep{depi}));
        if length(depthmatch) > 3
            % Try to recombine the files, starting with 2
            files2combo = mlist(depthmatch(2:end));
            combineAOfiles(matdest, files2combo)
        end
    end
end



% Generate the "Key" file that's needed for the rest of the processing
% Requires AO files are downloaded and converted, 

function AutoGenKey(subjstr)
    disp('Generating File Key');
    dddir = getenv('DBSDATADIR');

    % Get the AO task files
    AOmatdir = [dddir '\' subjstr '\AO_Case_Data\matfiles'];
    AOmatfiles = dir([AOmatdir '/*.mat']);
    AOmatfiles = {AOmatfiles.name};
    Depths = [];
    Fs = [];
    di = 1;
    % Scan the alpha omega mat files
    for i = 1:length(AOmatfiles)
        % Load each file
        AO = load([AOmatdir '/' AOmatfiles{i}]);
        % Check if digital in 1 and 2 exist and contain data and recording is
        % longer than 2 minutes
        if isfield(AO, {'CDIG_IN_1_Up', 'CDIG_IN_2_Up'})
            if ~isempty(AO.CDIG_IN_1_Up) && ~isempty(AO.CDIG_IN_2_Up) && (AO.CDIG_IN_1_TimeEnd - AO.CDIG_IN_1_TimeBegin) > 90
                Depths(di) = str2num(cell2mat(extractBetween(AOmatfiles{i}, 'D','F')));
                Fs(di) = str2num(cell2mat(extractBetween(AOmatfiles{i}, 'F', '.mat')));
                di = di+1;
    %             disp(AOmatfiles{i});
            end
        elseif isfield(AO, {'CDIG_IN_2_Up', 'CDIG_IN_3_Up'}) % Also allow it to be in digital inputs 2 and 3 (whoops)
            if ~isempty(AO.CDIG_IN_2_Up) && ~isempty(AO.CDIG_IN_3_Up) && (AO.CDIG_IN_2_TimeEnd - AO.CDIG_IN_2_TimeBegin) > 90
                Depths(di) = str2num(cell2mat(extractBetween(AOmatfiles{i}, 'D','F')));
                Fs(di) = str2num(cell2mat(extractBetween(AOmatfiles{i}, 'F', '.mat')));
                di = di+1;
    %             disp(AOmatfiles{i});
            end
        end
    end
    nAO = length(Depths);
    Depths = fliplr(Depths);
    Fs = fliplr(Fs);

    % Get task session numbers
    taskfolds = dir([dddir '\' subjstr '\' subjstr '_Task\Session*']);
    nTask = length(taskfolds);
    Tasks = nan(nTask,1);
    for i = 1:nTask
        sessparts = strsplit(taskfolds(i).name, ' ');
        Tasks(i) = str2num(sessparts{2});
    end

    % Get video session numbers
    vidfolds = dir([dddir '\' subjstr '\' subjstr '_Compressed_Video\' subjstr '\session*']);
    nVid = length(vidfolds);
    Vids = nan(nVid,1);
    for i = 1:nVid
        Vids(i) = str2num(vidfolds(i).name(end-2:end));
    end

    if (nAO ~= nTask) || (nTask ~= nVid) || (nAO ~= nVid)
        disp('UNMATCHED NUMBERS OF AO FILES / TASK / VIDEO, MAKE KEYFILE BY HAND');
        fprintf('Identified Depths and Fs:\n')
        for i = 1:nAO
            fprintf([num2str(Depths(i)) ' ' num2str(Fs(i)) '\n']);
        end
        fprintf('Identified Task Sessions:\n')
        for i = 1:nTask
            fprintf([num2str(Tasks(i)) '\n']);
        end
        fprintf('Identified Video Sessions:\n');
        for i = 1:nVid
            fprintf([num2str(Vids(i)) '\n']);
        end
        fid = fopen([dddir '\' subjstr '\' subjstr '_Key'], 'wt');
        fprintf(fid, 'AODepth\tF\tTas\tVid\n');
        fclose(fid);
        
        input('Go make that key by hand, then enter ''k'' here when done: ');
    else
        % Make Key
        fid = fopen([dddir '\' subjstr '\' subjstr '_Key'], 'wt');
        fprintf(fid, 'AODepth\tF\tTas\tVid\n');
        for i = 1:nAO
            fprintf(fid, [num2str(Depths(i)) '\t' num2str(Fs(i)) '\t' num2str(Tasks(i)) '\t' num2str(Vids(i))]);
            if i < nAO
                fprintf(fid, '\n');
            end
        end
        fclose(fid);
        disp(['Key File Generated for ' subjstr]);
    end
end


% Convert AO files to .plx files!
% Using write_plx from the plexon SDK.
% write_plx only writes waveform snippets. Therefore we must threshold here
% in this code and get snippets and threshold crossing times, then export
% those.

function Batch_AO2plx(subjstr)

    dddir = getenv('DBSDATADIR');

    AOmatdir = [dddir '\' subjstr '\AO_Case_Data\matfiles'];

    % Load key data
    keyfilename = [dddir '\' subjstr '\' subjstr '_Key'];
    K = tdfread(keyfilename);

    nfi = length(K.AODepth);

    % Discover AO mat file naming format
    AOmatfiles = dir([AOmatdir '/*.mat']);
    AO1 = AOmatfiles(1).name;
    AOprefix = strsplit(AO1,'D');
    AOprefix = AOprefix{1};

    savedir = [getenv('PLXDIR') '\' subjstr];
    if ~exist(savedir)
        mkdir(savedir)
    end

    for fii = 1:nfi
        % Get the AO mat data
        AOfile = dir([AOmatdir '\*D' num2str(K.AODepth(fii), '%0.3f') 'F' num2str(K.F(fii), '%0.4i') '.mat']);
        AOfile = [AOmatdir '\' AOfile.name];

        disp(['Doing ' AOfile]);
        AOmat_to_plx( AOfile, savedir);
    end
end

% Convert AO .mat files to .plx files!
%
% Takes the CSPK fields from an AO .mat file, finds threshold crossings,
% extracts snippets and writes them to a .plx file. Note that write_plx.m
% only lets you write one channel, so if the AO file has >1 CSPK fields,
% they will be written to separate .plx files.
%
% Required software: Plexon OmniPlex and MAP Offline MATLAB SDK (can be
% found at https://plexon.com/software-downloads/#software-downloads-SDKs).
% 
% After downloading the SDK, add to your MATLAB path.
%
% Required inputs:
% - AOfile: location of the Alpha Omega .mat file you wish to convert (must
%   have already been processed with AO Converter
% - savedir: location of folder to save the output .plx file
%
% Optional inputs:
% - thresh: Threshold in multiples of std. dev. to trigger snippet
%   extraction (can be positive or negative). Default: -2.5
% - snipt: The length of each snippet in ms. Default: 0.8
% - pbefore: The fraction of the snippet to occur before the threshold
% crossing. Default: 0.25
% - adaptive: Boolean flag determining if adaptive snippets should be used.
% Default: 0
% - adaptive_t: Time window in seconds to calculate the threshold over.
% Default: 10
% - adaptive_step: Step size in seconds to recalculate the threshold.
% Default: 2
% - align_peaks: Boolean flag determining if snippets will be recentered to
% align peaks. Can be 'none', 'pos' (align peaks), or 'neg' (align troughs).
% Defualt: 'none'
% - rescale: Boolean flag determining if waveforms will be rescaled,
% perhaps helping with the whole 

% Rex Tien 6/29/21
% 
function AOmat_to_plx( AOfile, savedir, varargin)

    plxsdkdir = getenv('PLXSDKDIR');
    addpath(genpath(plxsdkdir));
    
    thresh_def = -3;
    snipt_def = 0.8;
    pbefore_def = 0.25;
    adaptive_def = 0;
    adaptive_t_def = 10;
    adaptive_step_def = 2;
    align_peaks_def = 'none';
    p = inputParser;
    isfi = @(x) ~isempty(x) && isfile(x);
    isdi = @(x) ~isempty(x) && exist(x, 'dir');
    addRequired(p, 'AOfile', isfi);
    addRequired(p, 'savedir', isdi);
    addOptional(p, 'thresh', thresh_def);
    addOptional(p, 'snipt', snipt_def);
    addOptional(p, 'pbefore', pbefore_def);
    addOptional(p, 'adaptive', adaptive_def);
    addOptional(p, 'adaptive_t', adaptive_t_def);
    addOptional(p, 'adaptive_step', adaptive_step_def);
    addOptional(p, 'align_peaks', align_peaks_def);
    parse(p, AOfile, savedir, varargin{:});
    
    invar = fieldnames(p.Results);
    for vari = 1:length(invar)
        eval([invar{vari} ' = p.Results.' invar{vari} ';']);
    end
    
    if ~any(strcmp(align_peaks, {'none', 'pos', 'neg'}))
        disp('align_peaks must have a value of none, pos or neg. Exiting');
        return;
    end
    
    AOfile = p.Results.AOfile;
    savedir = p.Results.savedir;
    
    A = load(AOfile);
    
    AOfileonly = strsplit(AOfile, '\');
    AOfileonly = AOfileonly{end};
    
    if isfield(A,'CSPK_01_KHz')
        freq = A.CSPK_01_KHz*1000;
    elseif isfield(A,'CSPK_02_KHz')
        freq = A.CSPK_02_KHz*1000;
    elseif isfield(A,'CSPK_03_KHz')
        freq = A.CSPK_03_KHz*1000;
    end
    
    npw = ceil((snipt/1000)*freq);
    snipt = npw*(1/freq); % Actual snippet length
    nbefore = round(pbefore*npw); % number of points before crossing
    nafter = npw-nbefore;
    
    % Do the thresholding and extract snippets
    for i = 1:5
        if isfield(A, ['CSPK_' num2str(i,'%0.2i')])
            
            savefile = [savedir '\' AOfileonly(1:end-4) '_' num2str(i) '.plx'];
            
            spkstr = ['CSPK_' num2str(i, '%0.2i')];
            spkv = double(A.(spkstr))*A.([spkstr '_BitResolution'])/A.([spkstr '_Gain']);
            allts = (0:(length(spkv)-1))*(1/freq);
            
            if ~adaptive % Do not adaptive?
                threshv = thresh*std(spkv);

                % Find threshold crossing times
                if sign(thresh) == -1
                    ts = allts(find(spkv(1:end-1) > threshv & spkv(2:end) < threshv));
                elseif sign(thresh) == 1
                    ts = allts(find(spkv(1:end-1) < threshv & spkv(2:end) > threshv));
                end
            else % Do adaptive?
                step_starts = 0:adaptive_step:allts(end);
                step_mids = step_starts + adaptive_step/2;
                nstep = length(step_starts);
                ts = [];
                for stepi = 1:nstep
                    threshv = thresh*std(spkv( allts>(step_mids(stepi)-adaptive_t/2) & allts<(step_mids(stepi)+adaptive_t/2) ));
                    % Find threshold crossing times
                    spkvwin = spkv(allts>(step_mids(stepi)-adaptive_step/2) & allts<(step_mids(stepi)+adaptive_step/2));
                    tswin = allts(allts>(step_mids(stepi)-adaptive_step/2) & allts<(step_mids(stepi)+adaptive_step/2));
                    if sign(thresh) == -1
                        ts = [ts, tswin(find(spkvwin(1:end-1) > threshv & spkvwin(2:end) < threshv))];
                    elseif sign(thresh) == 1
                        ts = [ts, tswin(find(spkvwin(1:end-1) < threshv & spkvwin(2:end) > threshv))];
                    end
                end
            end            
            
            % Eliminate times too close to the start or end
            while ts(1) <= nbefore*(1/freq)
                ts(1) = [];
            end
            while ts(end) >= allts(end)-(nafter+1)*(1/freq)
                ts(end) = [];
            end
            
            % Eliminate times within 1 snipt of the previous crossing
            tsi = 1;
            while tsi <= length(ts)-1
                if ts(tsi+1)-ts(tsi) <= snipt
                    ts(tsi) = [];
                else
                    tsi = tsi+1;
                end
            end
            
            % Get the snippets
            n = length(ts);
            wave = nan(n,npw);
            if mod(npw,2) == 0
                peakbefore = npw/2-1;
                peakafter = npw/2;
            else
                peakbefore = floor(npw/2);
                peakafter = floor(nwp/2);
            end
            
            for j = 1:n
                tpos = find(allts == ts(j));
                wave(j,:) = spkv(tpos-nbefore+1:tpos+nafter);
                if strcmp(align_peaks, 'pos')
                    [~,peaki] = max(wave(j,:));
                    ppos = tpos+peaki-nbefore+1;
                    wave(j,:) = spkv(ppos-peakbefore:ppos+peakafter);
                elseif strcmp(align_peaks, 'neg')
                    [~,peaki] = min(wave(j,:));
                    ppos = tpos+peaki-nbefore+1;
                    wave(j,:) = spkv(ppos-peakbefore:ppos+peakafter);
                end
            end
            
            
            % Write the plx file
            eval(['write_plx(savefile, ' num2str(i) ', freq, npw, n, ts, wave, zeros(n,1))']);
            
        end
    end
end

% Auto-generate DLC python scripts, run first couple steps

function GenRunDLCPy(subjstr, ndim)
    % Where is the Key file?
    dddir = getenv('DBSDATADIR');
    dlcdir = getenv('DLCDIR');
    pydir = getenv('ANACONDAPYDIR');

    Keyfilename = [dddir '/' subjstr '/' subjstr '_Key'];
    Key = tdfread(Keyfilename);    
    
    rawviddir = [dddir '\' subjstr '\' subjstr '_Compressed_Video\' subjstr '\'];
    rawviddirb = strrep(rawviddir(1:end-1), '\', '/');

    nsess = length(Key.Vid);
    
    for si = 1:nsess
        sesses(si).name = ['session' num2str(Key.Vid(si),'%0.3i')];
    end
    s1dir = dir([rawviddir sesses(1).name '\2*']);

    viddate = s1dir(1).name(1:8);

    DLCdir = [dlcdir '\' subjstr '_' num2str(ndim) 'D\'];
    if ~exist(DLCdir)
        mkdir(DLCdir);
    end
    cd(DLCdir)

    if ndim == 3
        % possible camera names
        camnames = {'left', 'top', 'right'};
        ncam = length(camnames);
    end

    vidstr = ['['];
    for sessi = 1:nsess
        if ndim == 2
            vidstr = [vidstr 'viddir+''/' sesses(sessi).name '/' viddate '_' subjstr '_' sesses(sessi).name '_topCam-0000.mp4'''];
            if sessi == nsess
                vidstr = [vidstr ']'];
            else
                vidstr = [vidstr ', '];
            end
        elseif ndim == 3
            vidstub = ['viddir+''/' sesses(sessi).name '/' viddate '_' subjstr '_' sesses(sessi).name];
            for cami = 1:ncam
                vidstr = [vidstr vidstub '_' camnames{cami} 'Cam-0000.mp4'''];
                if (sessi == nsess) && (cami == ncam)
                    vidstr = [vidstr ']'];
                else
                    vidstr = [vidstr ', '];
                end
            end
        end
    end

    % Make and run init.py
    fid = fopen([subjstr '_init.py'], 'wt');
    fprintf(fid, [...
        'import deeplabcut\n' ...
        'task = ''' subjstr '''\n' ...
        'experimenter = ' '''NERDCO''' '\n' ...
        'viddir = ''' rawviddirb '''\n' ...
        'video = ' vidstr '\n' ...
        'path_config_file = deeplabcut.create_new_project(task,experimenter,video,copy_videos=True)']);
    fclose(fid);

    eval(['!' pydir ' ' subjstr '_init.py']);

    % Find and edit the conf file
    if ndim == 2
        projdir = dir([subjstr '-NERDCO*']);
        confdir = [DLCdir projdir(1).name '\config.yaml'];
        conf = readlines(confdir);
        bpline = find(strcmp('bodyparts:', conf));
        conf(bpline+1) = '- fingerPoint';
        conf(bpline+2:bpline+4) = [];
        nfpline = find(strcmp('numframes2pick: 20', conf));
        conf(nfpline) = ['numframes2pick: ' num2str(ceil(200/nsess))];
        skline = find(strcmp('skeleton:', conf));
        conf(skline) = 'skeleton: ';
        conf(skline+1:skline+4) = [];
        dlline = find(strcmp('dotsize: 12', conf));
        conf(dlline) = 'dotsize: 3';
        fid = fopen(confdir, 'w');
        fwrite(fid, strjoin(conf, '\n'));
        fclose(fid);
    elseif ndim == 3
        projdir = dir([subjstr '-NERDCO*']);
        confdir = [DLCdir projdir(1).name '\config.yaml'];
        conf = readlines(confdir);
        
        % Pad out conf
        for padi = 1:20
            conf(end+1) = '';
        end
    
        % Fix bodyparts
        start0line = find(strcmp('start: 0', conf));
        starton = conf(start0line:end);
        bpline = find(strcmp('bodyparts:', conf));
        conf(bpline+1) = '- fingerTip2';
        conf(bpline+2) = '- MCP2';
        conf(bpline+3) = '- MCP5';
        conf(bpline+4) = '- wristCenter';
        conf(bpline+5) = '- stablePoint';
        conf(bpline+6:end) = starton(1:end-1);
        
        % Fix NFrames2Pick
        nfpline = find(strcmp('numframes2pick: 20', conf));
        conf(nfpline) = ['numframes2pick: ' num2str(50)];
        
        % Fix Skeleton
        scolline = find(strcmp('skeleton_color: black', conf));
        scollon = conf(scolline:end);
        skline = find(strcmp('skeleton:', conf));
        conf(skline+1) = '- - fingerTip2';
        conf(skline+2) = '  - MCP2';
        conf(skline+3) = '- - fingerTip2';
        conf(skline+4) = '  - MCP5';
        conf(skline+5) = '- - fingerTip2';
        conf(skline+6) = '  - wristCenter';
        conf(skline+7) = '- - fingerTip2';
        conf(skline+8) = '  - stablePoint';
        conf(skline+9) = '- - MCP2';
        conf(skline+10) = '  - MCP5';
        conf(skline+11) = '- - MCP2';
        conf(skline+12) = '  - wristCenter';
        conf(skline+13) = '- - MCP2';
        conf(skline+14) = '  - stablePoint';
        conf(skline+15) = '- - MCP5';
        conf(skline+16) = '  - wristCenter';
        conf(skline+17) = '- - MCP5';
        conf(skline+18) = '  - stablePoint';
        conf(skline+19) = '- - wristCenter';
        conf(skline+20) = '  - stablePoint';
        conf(skline+21:end) = scollon(1:end-16);
        
        % Fix DotSize
        dlline = find(strcmp('dotsize: 12', conf));
        conf(dlline) = 'dotsize: 3';
        
        % Trim extra
        lastline = find(diff(strcmp('', conf))==1,1,'last');
        conf = conf(1:lastline+1);
        
        fid = fopen(confdir, 'w');
        fwrite(fid, strjoin(conf, '\n'));
        fclose(fid);
    end

    confdirb = strrep(confdir, '\', '/');

    newviddir = [DLCdir projdir(1).name '\videos'];
    newviddirb = strrep(newviddir, '\', '/');

    newvidstr = ['['];
    for sessi = 1:nsess
        if ndim == 2
            newvidstr = [newvidstr 'viddir+''/' viddate '_' subjstr '_' sesses(sessi).name '_topCam-0000.mp4'''];
            if sessi == nsess
                newvidstr = [newvidstr ']'];
            else
                newvidstr = [newvidstr ', '];
            end
        elseif ndim == 3
            newvidstub = ['viddir+''/' viddate '_' subjstr '_' sesses(sessi).name];
            for cami = 1:ncam
                newvidstr = [newvidstr newvidstub '_' camnames{cami} 'Cam-0000.mp4'''];
                if (sessi == nsess) && (cami == ncam)
                    newvidstr = [newvidstr ']'];
                else
                    newvidstr = [newvidstr ', '];
                end
            end
        end
    end

    % Make prep_labels.py
    fid = fopen([subjstr '_prep_labels.py'], 'wt');
    fprintf(fid, [...
        'import deeplabcut\n' ...
        'path_config_file = ''' confdirb '''\n' ...
        'deeplabcut.extract_frames(path_config_file, userfeedback=False)\n']);
    fclose(fid);

    % Make label.py
    fid = fopen([subjstr '_label.py'], 'wt');
    fprintf(fid, [...
        'import deeplabcut\n' ...
        'path_config_file = ''' confdirb '''\n' ...
        'deeplabcut.label_frames(path_config_file)\n' ...
        'deeplabcut.check_labels(path_config_file)\n' ...
        'deeplabcut.create_training_dataset(path_config_file)\n']);
    fclose(fid);
    
    % Make train.py
    fid = fopen([subjstr '_train.py'], 'wt');
    fprintf(fid, [...
        'import deeplabcut\n' ...
        'path_config_file = ''' confdirb '''\n' ...
        'deeplabcut.train_network(path_config_file, allow_growth=True)\n']);
    fclose(fid);
    
    % Make eval-analyze.py
    fid = fopen([subjstr '_eval-analyze.py'], 'wt');
    fprintf(fid, [...
        'import deeplabcut\n' ...
        'path_config_file = ''' confdirb '''\n' ...
        'if __name__ == ''__main__'':\n' ...
        '    deeplabcut.evaluate_network(path_config_file)\n' ...
        '    viddir = ''' newviddirb '''\n' ...
        '    videofile_path = ' newvidstr '\n' ...
        '    deeplabcut.analyze_videos(path_config_file, videofile_path, videotype=''.mp4'', save_as_csv=True)\n' ...
        '    deeplabcut.create_labeled_video(path_config_file, videofile_path, videotype=''.mp4'')\n' ...
        '    deeplabcut.plot_trajectories(path_config_file, videofile_path)']);
    fclose(fid);
    
    % Run prep_labels.py
    eval(['!' pydir ' ' subjstr '_prep_labels.py']);

    disp('OK, please go run _label.py, _train.py and _eval-analyze.py from anaconda prompt. Sorry I can''t do it from here');
end

function combineAOfiles(destdir, f2c)
    % Check that first file is good, then start combining til the end of
    % the list of files to combine
    if ~isempty(f2c)
        load([destdir '\' f2c{1}]);
        if exist('CDIG_IN_1_Up','var') && exist('CDIG_IN_2_Up','var')
            if ~isempty(CDIG_IN_1_Up) && ~isempty(CDIG_IN_2_Up) && (CDIG_IN_1_TimeEnd - CDIG_IN_1_TimeBegin) > 30
                for combi = 2:length(f2c)
                    B = load([destdir '\' f2c{combi}]);
                    if isfield(B, {'CDIG_IN_1_Up', 'CDIG_IN_2_Up'})
                        for digi = 1:5
                            if exist(['CDIG_IN_' num2str(digi) '_Up'],'var')
                                eval(['CDIG_IN_' num2str(digi) '_Up = [CDIG_IN_' num2str(digi) '_Up, B.CDIG_IN_' num2str(digi) '_Up + (B.CDIG_IN_' num2str(digi) '_TimeBegin-CDIG_IN_' num2str(digi) '_TimeBegin)*(1000*CDIG_IN_' num2str(digi) '_KHz)];'])
                                eval(['CDIG_IN_' num2str(digi) '_Down = [CDIG_IN_' num2str(digi) '_Down, B.CDIG_IN_' num2str(digi) '_Down + (B.CDIG_IN_' num2str(digi) '_TimeBegin-CDIG_IN_' num2str(digi) '_TimeBegin)*(1000*CDIG_IN_' num2str(digi) '_KHz)];'])
                                eval(['CDIG_IN_' num2str(digi) '_TimeEnd = B.CDIG_IN_' num2str(digi) '_TimeEnd;']);
                            end
                        end
                        
                        for reci = 1:5
                            if exist(['CLFP_' num2str(reci,'%0.2i')], 'var')
                                eval(['CLFP_' num2str(reci,'%0.2i') ' = [CLFP_' num2str(reci,'%0.2i') ', B.CLFP_' num2str(reci,'%0.2i') '];']);
                                eval(['CLFP_' num2str(reci,'%0.2i') '_TimeEnd = B.CLFP_' num2str(reci,'%0.2i') '_TimeEnd;']);
                                
                                eval(['CMacro_LFP_' num2str(reci,'%0.2i') ' = [CMacro_LFP_' num2str(reci,'%0.2i') ', B.CMacro_LFP_' num2str(reci,'%0.2i') '];']);
                                eval(['CMacro_LFP_' num2str(reci,'%0.2i') '_TimeEnd = B.CMacro_LFP_' num2str(reci,'%0.2i') '_TimeEnd;']);
                                
                                eval(['CMacro_RAW_' num2str(reci,'%0.2i') ' = [CMacro_RAW_' num2str(reci,'%0.2i') ', B.CMacro_RAW_' num2str(reci,'%0.2i') '];']);
                                eval(['CMacro_RAW_' num2str(reci,'%0.2i') '_TimeEnd = B.CMacro_RAW_' num2str(reci,'%0.2i') '_TimeEnd;']);
                                
                                eval(['CRAW_' num2str(reci,'%0.2i') ' = [CRAW_' num2str(reci,'%0.2i') ', B.CRAW_' num2str(reci,'%0.2i') '];']);
                                eval(['CRAW_' num2str(reci,'%0.2i') '_TimeEnd = B.CRAW_' num2str(reci,'%0.2i') '_TimeEnd;']);
                                
                                eval(['CSPK_' num2str(reci,'%0.2i') ' = [CSPK_' num2str(reci,'%0.2i') ', B.CSPK_' num2str(reci,'%0.2i') '];']);
                                eval(['CSPK_' num2str(reci,'%0.2i') '_TimeEnd = B.CSPK_' num2str(reci,'%0.2i') '_TimeEnd;']);
                            end
                        end
                    end
                    % Get rid of file after combining
                    delete([destdir '\' f2c{combi}]);
                end
                save([destdir '\' f2c{1}],'-regexp', '^(?!(destdir|f2c|B|combi|digi|reci)$).');
            else
                disp([f2c{1} ' is no good, iterating...']);
                combineAOfiles(destdir, f2c(2:end));
            end
        else
            disp([f2c{1} ' is no good, iterating...']);
            combineAOfiles(destdir, f2c(2:end));
        end
    end
end