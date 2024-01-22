% Function to perform fixing for skipped frames in the DLC video recording.
%
% Rex Tien 12/14/23

% Frame skip fixing:
%
% If the time vector looks fine, don't alter it.
%
% 1. Detect when the difference between frame times is a multiple of the
% framerate
%
% 2. Linear interpolation between frame times on either side of skips to
% fill in the time vector so that in the end, times occur at approximately
% the correct framerate.
%
% 3. Insert NaNs into the position matrix to represent missed tracking
% during the skipped frames

% Required inputs:
% gapktime: original time vector from DLC alignment (in seconds)
% rawpos: original positions matrix from DLC output (NxD matrix. As many
%       positions can be sent in as desired, each dimension as a column)
% framerate: the intended framerate of the video recording

% Optional settings:
% verbose: bool - Print things along the way? (default = true)

% Outputs:
% ktime: Time vector expanded to accomodate skipped frames
% pos: Position matrix with nans inserted where skipped frames existed

function [ktime, pos] = fix_frameskips(gapktime, rawpos, framerate, varargin)
    
    % Default values:
    V.verbose = true;
    
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
    
    boolchecks = {'verbose'};
    for booli = 1:length(boolchecks)
        if ~isa(V.(boolchecks{booli}), 'logical')
            error([boolchecks{booli} ' must be a boolean value (true or false)']);
        end
    end
    
    if size(gapktime,1) ~= size(rawpos,1)
        error('input time and raw position must have the same number of samples (rows)');
    end
    
    if any(isnan(gapktime))
        error('NaN value(s) detected in input time')
    end
    
    ndim = size(rawpos,2);
    if ndim<=0
        error('position matrix appears to be empty');
    end
    
    framestep = 1/framerate;
    
    % Do the frame skip fix if necessary
    dkt = diff(gapktime);
    % check for weird short gaps
    if any(unique(dkt) < 0.9*framestep)
        error('detected abnormally short time between frames - check time vector');
    end
    % Detect abnormally long gaps that would suggest missing video frames
    if any(unique(dkt) > 1.1*framestep) || any(unique(dkt) < 0.9*framestep)
        if V.verbose
            disp('Detected skipped frame(s), doing fix.');
        end
        % Convert time to samples by rounding
        fakeusetime = [0; cumsum(round(dkt/framestep))];
        fakestep = 1;
        dut = diff(fakeusetime);
    
        % Try to calculate the number of samples there should be
        nsampfix = round((fakeusetime(end)-fakeusetime(1))/fakestep)+1;
        fixusetime = nan(nsampfix,1);
        ktime = nan(nsampfix,1);
        addedpt = zeros(nsampfix,1);
        ndim = size(rawpos,2);
        pos = nan(nsampfix,ndim);
        fixcount = 1;
        for fixi = 1:length(dut)
            % if it's a good sample, store it and increment; if
            % it's bad store the current sample and the skipped
            % ones (not the next sample)
            nskipped = round(dut(fixi)/fakestep);
            % interpolate time during skipped frames
            addtime = interp1([fixi, fixi+nskipped], [fakeusetime(fixi), fakeusetime(fixi+1)], fixi:fixi+nskipped);
            fixusetime(fixcount:fixcount+nskipped-1) = addtime(1:end-1);
            addtimeorig = interp1([fixi, fixi+nskipped], [gapktime(fixi), gapktime(fixi+1)], fixi:fixi+nskipped);
            ktime(fixcount:fixcount+nskipped-1) = addtimeorig(1:end-1);
            % fill in nans in kinematics
            pos(fixcount:fixcount+nskipped-1,:) = [rawpos(fixi,:); nan(nskipped-1,ndim)];
            if (nskipped>1)
                addedpt(fixcount+1:fixcount+nskipped-1) = ones(nskipped-1,1);
            end
            fixcount = fixcount+nskipped;
        end
        fixusetime(end) = fakeusetime(end);
        ktime(end) = gapktime(end);
        pos(end,:) = rawpos(end,:);
        if V.verbose
            disp(['Fixed ' num2str(sum(addedpt)) ' missing frames']);
        end
    else
        % If no gaps, let it through
        if V.verbose
            disp('No skipped frames detected, returning');
        end
        ktime = gapktime;
        pos = rawpos;
    end