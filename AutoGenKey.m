% Generate the "Key" file that's needed for the rest of the processing
% Requires AO files are downloaded and converted, 

% subjstr = 'PD07';
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