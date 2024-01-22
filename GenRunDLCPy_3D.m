% Auto-generate DLC python scripts

% function GenRunDLCPy(subjstr)
    subjstr = 'VM31';

    % Where is the Key file?
    dddir = getenv('DBSDATADIR');
%     Keyfilename = [dddir '/' subjstr '/' subjstr '_Key'];
%     Key = tdfread(Keyfilename);
    Key.Vid = [1 2];
    
    
    rawviddir = [dddir '\' subjstr '\' subjstr '_Compressed_Video\' subjstr '\'];
    rawviddirb = strrep(rawviddir(1:end-1), '\', '/');
    pydir = 'C:\Users\RexBoxOne\anaconda3\envs\DLC-GPU\python.exe';

    nsess = length(Key.Vid);
    
    for si = 1:nsess
        sesses(si).name = ['session' num2str(Key.Vid(si),'%0.3i')];
    end
    s1dir = dir([rawviddir sesses(1).name '\2*']);

    viddate = s1dir(1).name(1:8);

    DLCdir = ['C:\Users\RexBoxOne\Documents\DeepLabCut\dbsdlc\' subjstr '_3D\'];
    if ~exist(DLCdir)
        mkdir(DLCdir);
    end
    cd(DLCdir)

    % possible camera names
    camnames = {'left', 'top', 'right'};
    ncam = length(camnames);
    
    vidstr = ['['];
    for sessi = 1:nsess
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

    % Make init.py
    fid = fopen([subjstr '_init.py'], 'wt');
    fprintf(fid, [...
        'import deeplabcut\n' ...
        'task = ''' subjstr '''\n' ...
        'experimenter = ' '''Rex''' '\n' ...
        'viddir = ''' rawviddirb '''\n' ...
        'video = ' vidstr '\n' ...
        'path_config_file = deeplabcut.create_new_project(task,experimenter,video,copy_videos=True)']);
    fclose(fid);

    eval(['!' pydir ' ' subjstr '_init.py']);

    % Find and edit the conf file
    projdir = dir([subjstr '-Rex*']);
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

    confdirb = strrep(confdir, '\', '/');

    newviddir = [DLCdir projdir(1).name '\videos'];
    newviddirb = strrep(newviddir, '\', '/');

    newvidstr = ['['];
    for sessi = 1:nsess
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