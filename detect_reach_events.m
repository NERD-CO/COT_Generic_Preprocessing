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
        if ~strcmp(vidnotes{1,4}, 'g')
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
