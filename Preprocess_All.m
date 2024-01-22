% Run all analysis steps!

clear all;
close all;

% Suppress the spline nan warning.
warning('off','MATLAB:chckxy:IgnoreNaN');

% Get important directories
datadir = getenv('PDGENERIC2DDATA');
dddir = getenv('DBSDATADIR');

% Set and retrieve project settings
SetSS_COT_STN(datadir);
Settings = load([datadir 'SS.mat']);

ndim = 2;

% regen = 'none';
regen = 'rkn';

for subi = 1:length(Settings.SS)
    SSi = Settings.SS{subi};

    % Find the data tartares
    tardir = [dddir SSi '\Data_Tartare_Full_' num2str(ndim) 'D\'];

    if ~isfolder(tardir)
        mkdir(tardir);
    end

    tarfiles = dir([tardir '*' SSi '*Full_' num2str(ndim) 'D.mat']);

    % Make the data tartares if they don't yet exist
    if isempty(tarfiles) || any(strcmp(regen,'tar'))
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
            % LATER YOU CAN RERUN IT WITH 0 ONCE YOU'VE SCREENED EVERYTHING
            [Reach] = detect_reach_events(SSi, fii, Time.ktime, Kin.SmoothG50.Vel, Tar.T, Settings.vidframerate, [annodir annofile(1).name], [dddir 'COT_ReachFind_Fixes.xlsx'], 0);
    
            % Save current stuff as a reach_kin_neu file
            disp('');
            N = Tar.N;
            N.SPK
%             % remove SPK to save space because I'm not using it here
%             N = rmfield(N,'SPK');
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