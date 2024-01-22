% Pull files from synology, make nice folder
% to do: PD06-10 VM10-14
% subjstr = 'PD07';
function PullFromSynology(subjstr)
    disp(['Pulling files for ' subjstr]);

%     syndir = 'S:\CenterOutTask';
    syndir = '\\som-nsg-rsrch5.ucdenver.pvt\documents\CenterOutTask';

    dddir = getenv('DBSDATADIR');

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
    cd 'C:\Program Files (x86)\AlphaOmega\Converter';
    for i = 1:length(mpxlist)
        fprintf([num2str(i,'%0.3i') '/' num2str(length(mpxlist),'%0.3i')]);
        mpxsc = [dddir '\' subjstr '\AO_Case_Data\' mpxlist(i).name];
        command=['Mapfile.exe =' mpxsc '=Matlab=' matdest '='];
        [status,cmdout] = dos(command);
        fprintf('\b\b\b\b\b\b\b');
    end
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