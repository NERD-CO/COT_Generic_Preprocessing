% Auto-generate DLC python scripts

function GenRunDLCPy(subjstr)
%     subjstr = 'VM31';

    % Where is the Key file?
    dddir = getenv('DBSDATADIR');
    Keyfilename = [dddir '/' subjstr '/' subjstr '_Key'];
    Key = tdfread(Keyfilename);
%     Key.Vid = [1 2];
    
    
    rawviddir = [dddir '\' subjstr '\' subjstr '_Compressed_Video\' subjstr '\'];
    rawviddirb = strrep(rawviddir(1:end-1), '\', '/');
    pydir = 'C:\Users\RexBoxOne\anaconda3\envs\DLC-GPU\python.exe';

    nsess = length(Key.Vid);
    
    for si = 1:nsess
        sesses(si).name = ['session' num2str(Key.Vid(si),'%0.3i')];
    end
    s1dir = dir([rawviddir sesses(1).name '\2*']);

    viddate = s1dir(1).name(1:8);

    DLCdir = ['C:\Users\RexBoxOne\Documents\DeepLabCut\dbsdlc\' subjstr '\'];
    if ~exist(DLCdir)
        mkdir(DLCdir);
    end
    cd(DLCdir)


    vidstr = ['['];
    for sessi = 1:nsess
        vidstr = [vidstr 'viddir+''/' sesses(sessi).name '/' viddate '_' subjstr '_' sesses(sessi).name '_topCam-0000.mp4'''];
        if sessi == nsess
            vidstr = [vidstr ']'];
        else
            vidstr = [vidstr ', '];
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

    confdirb = strrep(confdir, '\', '/');

    newviddir = [DLCdir projdir(1).name '\videos'];
    newviddirb = strrep(newviddir, '\', '/');

    newvidstr = ['['];
    for sessi = 1:nsess
        newvidstr = [newvidstr 'viddir+''/' viddate '_' subjstr '_' sesses(sessi).name '_topCam-0000.mp4'''];
        if sessi == nsess
            newvidstr = [newvidstr ']'];
        else
            newvidstr = [newvidstr ', '];
        end
    end

    % Make label.py
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
    
    eval(['!' pydir ' ' subjstr '_prep_labels.py']);

    disp('OK, please go run _label.py, _train.py and _eval-analyze.py from anaconda prompt. Sorry I can''t do it from here');