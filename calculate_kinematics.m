% Function to get higher order kinematics (velocity, acceleration and
% associated terms) from smooth position input.
%
% Rex Tien 12/14/23

% Kinematic processing
%
% Differentiation is done using [diff] function
%
% Smoothing is done using GaussSmooth_Arbitrary and is done zero times. The 
% smoothed position is used as input to generate 'smooth' higher order kinematics, 
% and raw is used to generate raw. Higher order kinematics are resampled
% back to position time using [spline].

% Required inputs:
% ktime: Aligned timestamps from DLC after frame skip fixing
% Pos: struct with fields:
%       Raw: frame-fixed positions matrix, not smoothed
%       Smooth.GXX: frame-fixed positions matrix, smoothed with Gaussian
%           kernel of width XX milliseconds, can have any number of GXX fields

% Outputs:
% Time: Struct with fields
%   ktime: same as input ktime
%
% Kin: Struct with kinematics resampled to position times (ktime)
%   Raw: Non-smoothed - recommended not to use as higher order terms will be very jumpy
%       pos
%       vel
%       speed
%       acc
%       accmag: magnitude of acceleration
%       accsign: magnitude of acceleration multiplied by the dot product between velocity and acceleration (positive when increasing the current acceleration, negative when decreasing it)
%   SmoothGXX: Kinematics calculated from smooth position (Smoothed with Gauss kernel of width XX milliseconds) with no re-smoothing
%       (same subfields as Kin.Raw)

function [Time, Kin] = calculate_kinematics(ktime, Pos)

    % Do some input checking
    
    if size(ktime,1) ~= size(Pos.Raw,1)
        error('input time and position must have the same number of samples (rows)');
    end
    
    if any(isnan(ktime))
        error('NaN value(s) detected in input time')
    end
    
    ndim = size(Pos.Raw,2);
    if ndim<=0
        error('position matrix appears to be empty');
    end
    
    smoothtypes = fieldnames(Pos);
    
    % Setting up
    ndim = size(Pos.Raw,2);
    Time.ktime = ktime;
    vtime = mean([ktime(1:end-1), ktime(2:end)], 2);
    atime = mean([vtime(1:end-1), vtime(2:end)], 2);

    % Get velocity, acc, resample back to ktime, calculate other kinematics
    for sti = 1:length(smoothtypes)
        stype = smoothtypes{sti};

        Kin.(stype).Pos = Pos.(stype);

        tempvel = diff(Kin.(stype).Pos,1)./diff(ktime);
        tempacc = diff(tempvel,1)./diff(vtime);

        Kin.(stype).Vel = nan(length(ktime),ndim);
        for dimi = 1:ndim
            Kin.(stype).Vel(:,dimi) = spline(vtime, tempvel(:,dimi), ktime);
        end
        Kin.(stype).Vel(isnan(Kin.(stype).Pos(:,1)),:) = NaN;
        Kin.(stype).Speed = vecnorm(Kin.(stype).Vel')';
        Kin.(stype).Dirs = Kin.(stype).Vel./(repmat(Kin.(stype).Speed,[1,ndim]));

        Kin.(stype).Acc = nan(length(ktime),ndim);
        for dimi = 1:ndim
            Kin.(stype).Acc(:,dimi) = spline(atime, tempacc(:,dimi), ktime);
        end
        Kin.(stype).Acc(isnan(Kin.(stype).Pos(:,1)),:) = NaN;
        Kin.(stype).AMag = vecnorm(Kin.(stype).Acc')';

        avdot = sum((Kin.(stype).Vel./Kin.(stype).Speed).*(Kin.(stype).Acc./Kin.(stype).AMag),2);
        Kin.(stype).ASign = Kin.(stype).AMag.*avdot;
    end
end