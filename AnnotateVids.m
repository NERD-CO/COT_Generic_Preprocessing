% This script will add trial numbers and target numbers to videos for
% debugging and bad trial exclusion!

function AnnotateVids(SS, ndim)

dddir = getenv('DBSDATADIR');

for subi = 1:length(SS)
    SSi = SS{subi}
    
    % Where to save everything?
    savedir = [dddir '\' SSi '\Annotated_Video'];
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    else
        movefile(savedir, [dddir '\' SSi '\Orig_Annotated_Video']);
        mkdir(savedir);
    end

    % Where are the labeled videos to annotate? If not found, look for a
    % folder in Documents / dbsdlc and move it
    subjdir = dir([dddir '\' SSi '\' SSi '-*-20*']);
    if isempty(subjdir)
        disp('DLC folder not found, looking for it in DBS DLC folder');
        subjdir = dir([dlcdir '\' SSi '\' SSi '-*-20*']);
        sdname = subjdir(end).name;
        copyfile([dlcdir '\' SSi '\' sdname], [dddir '\' SSi '\' sdname]);
        disp('DLC folder found and copied!');
    else
        sdname = subjdir(1).name;
    end
    Labeleddir = [dddir '\' SSi '\' sdname];
    
    % Where is the Key file?
    Keyfilename = [dddir '\' SSi '\' SSi '_Key'];
    
    % Read in the Key file
    Key = tdfread(Keyfilename);
    
    nfiles = length(Key.AODepth);
    disp(['Found ' num2str(nfiles) ' total files']);
    
    for fii = 1:nfiles
        fprintf(['Doing session ' num2str(fii) ' (video session ' num2str(Key.Vid(fii)) ')\n']);
        
        labvid = dir([Labeleddir '\videos\*_' SSi '*_session' num2str(Key.Vid(fii),'%0.3i') '_topCam*labeled.mp4']);
        labviddir = [Labeleddir '\videos\' labvid(1).name];
        thissavevid = [savedir '\' labvid(1).name];

        v = VideoReader(labviddir);
        fprintf('Reading video frames... ');
        if (strcmp(SSi, 'VM27') && fii == 2) || (strcmp(SSi, 'VM29') && fii == 2) || (strcmp(SSi, 'VM29') && fii == 3)
            allframes = zeros([540,720,3,8500], 'uint8');
            ticker = 0;
            while hasFrame(v)
                if mod(ticker,4) == 0
                    allframes(:,:,:,(ticker/4)+1) = readFrame(v);
                    ticker = ticker+1;
                else
                    trashframe = readFrame(v);
                    ticker = ticker+1;
                end
            end
        else
            allframes = read(v);
        end
        fprintf('Done!\n');
        nframevid = v.NumFrames;

        [ktime, T] = Align_TTLs_And_Frames(SSi, fii, Key, ndim);

        % Downsample!
        if ~(strcmp(SSi, 'VM27') && fii==2) && ~(strcmp(SSi, 'VM29') && fii==2) && ~(strcmp(SSi, 'VM29') && fii==3)
            allframes = allframes(:,:,:,1:4:v.NumFrames);
        end
        dstime = ktime(1:4:nframevid);
        
        fprintf('Adding text... ');
        outtimes = T.Time.StartOfMovement;
        backtimes = T.Time.CenterIsShown;
        ntri = min([length(outtimes), length(backtimes)]);
        for ti = 1:ntri
            outframes = find(dstime >= outtimes(ti) & dstime < backtimes(ti));
            for fi = 1:length(outframes)
                allframes(:,:,:,outframes(fi)) = insertText(allframes(:,:,:,outframes(fi)),[100 100],['Reach ' num2str(ti) ', Targ ' num2str(T.Target_idx(ti)) ', Out'],'FontSize',24,'BoxColor','black','BoxOpacity',1,'TextColor','white');
            end
            if ti < ntri
                backframes = find(dstime >= backtimes(ti) & dstime < outtimes(ti+1));
            else
                backframes = find(dstime >= backtimes(ti));
            end
            for fi = 1:length(backframes)
                allframes(:,:,:,backframes(fi)) = insertText(allframes(:,:,:,backframes(fi)),[100 100],['Reach ' num2str(ti) ', Targ ' num2str(T.Target_idx(ti)) ', Back'],'FontSize',24,'BoxColor','black','BoxOpacity',1,'TextColor','white');
            end
            
        end
        fprintf('Done!\n');
        
        fprintf('Saving video... ');
        vw = VideoWriter([thissavevid(1:end-4) '_annotated'], 'MPEG-4');
        open(vw);
        writeVideo(vw,allframes);
        close(vw);
        fprintf('Done!\n');
        
        % Also create an excel file to note all the bad / questionable
        % trials
        sheet{1,1} = 'Trial';
        sheet{1,2} = 'Direction';
        sheet{1,3} = 'Judgement (blank = good, x = bad or no reach, l = late reach, q = questionable reach)';
        for j = 1:ntri
            sheet{1+(2*(j-1)+1),1} = j;
            sheet{1+(2*(j-1)+2),1} = j;
            sheet{1+(2*(j-1)+1),2} = 'Out';
            sheet{1+(2*(j-1)+2),2} = 'Back';
        end
        xlswrite([thissavevid(1:end-4) '_notes.xlsx'], sheet);
    end
end