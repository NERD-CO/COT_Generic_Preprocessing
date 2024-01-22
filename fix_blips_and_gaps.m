% Function to process the raw position data by fixing blips and gaps and
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