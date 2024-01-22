% Checks TTL and camera frame alignment!

function CheckAlignment(SS, ndim)
    
    dddir = getenv('DBSDATADIR');

    for subi = 1:length(SS)
        SSi = SS{subi}
        
        % Where is the Key file?
        Keyfilename = [dddir '\' SSi '\' SSi '_Key'];
        
        % Read in the Key file
        Key = tdfread(Keyfilename);
        
        nfiles = length(Key.AODepth);
        disp(['Found ' num2str(nfiles) ' total files']);
    
        for fii = 1:nfiles
            fii
   
            [ktime, ~, Warning, Plt] = Align_TTLs_And_Frames(SSi, fii, Key, ndim);

            % Plot
            figure(1);
    %         figure;
            subplot(2,2,1);
            plot(Plt.AO_Task_TTL_times, 'r', 'Marker', '.');
            hold on;
            plot(Plt.TaskConv,'b');
            
            subplot(2,2,2);
            plot(diff(Plt.AO_Task_TTL_times), 'r', 'Marker', '.');
            hold on;
            plot(diff(Plt.Task_TTL_times), 'b');
            
            % How many frames before aligni?
            if Plt.camAligni - (Plt.kinAligni-1) > 0
                camClip = Plt.camTime(Plt.camAligni-(Plt.kinAligni-1):end);
            else
                camClip = [nan(1,Plt.kinAligni-Plt.camAligni) Plt.camTime];
            end
                   
            figure(1);
            subplot(2,2,3);
            plot(camClip, 'r', 'Marker', '.');
            hold on;
            plot(ktime, 'b')
            
            figure(1);
            subplot(2,2,4);
            plot(diff(camClip), 'r', 'Marker', '.');
            hold on;
            plot(diff(ktime), 'b');
            
            sgtitle([SSi ' ' num2str(fii) '/' num2str(nfiles) Warning]);
    
            while ~waitforbuttonpress
            end
            
        end
    end
    close all;
end