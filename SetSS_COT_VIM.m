function SetSS(savedir)

% Here is where to set major analysis parameters and settings!

% SS is subject IDs
SS = {'VM02', 'VM03', 'VM05', 'VM06', 'VM07', 'VM08', 'VM09', 'VM10', 'VM12', 'VM14', 'VM16', 'VM17', 'VM18', 'VM19', 'VM20', 'VM21', 'VM22', 'VM23', 'VM25', 'VM26', 'VM27', 'VM28', 'VM29', 'VM30', 'VM31'};
% SC is subject codes - one per each unique person
SC = [1,        1,         2,       3,    3,      2,      4,      5,      6,      7,      8,      8,      9,      9,    10,     10,     11,     12,     13,     13,     14,     14,     15,     16,     16];
% NN is the public facing subject code for this project (S.##.Operation_Hemisphere)
NN = {'S.01.L', 'S.01.R',  'S.02.L','S.03.L', 'S.03.R', 'S.02.R', 'S.04.R', 'S.05.L', 'S.06.L', 'S.07.L', 'S.08.L', 'S.08.R', 'S.09.R', 'S.09.L', 'S.10.L', 'S.10.R', 'S.11.L', 'S.12.L', 'S.13.L', 'S.13.R', 'S.14.R', 'S.14.L', 'S.15.R', 'S.16.L', 'S.16.R'};

kinname = repmat({'fingerTip2'}, [length(SS),1]);
kinname{strcmp(SS,'VM05')} = 'fingerTip';

ismale = [1 1 1 0 1 0 0 1 1 1 0 0 1 1 1 1];
ages = [71 71 74 71 75 50 73 68 62 72 65 68 73 80 77 73];

% Likelihood cutoff for rejecting DLC tracking
likecutoff = 0.95;

% Also set the video framerate and framestep for each subject
vidframerate = 120; % In fps

% Settings for screening out rateskin files / neurons
reachmin = 25;
FRmin = 1;

% Set the sliding renormalization parameters
slidewin = 45; % Sliding window in seconds
slidepad = (slidewin/2)./(1/vidframerate); % Amount to pad for sliding

% Significance vals to use for... everything
alph2do = [0.05 0.01];
cutoff = 5;

% Padding for rates before and after first and last reaches
msratepad = 10;

% Raster time window size
Raster.winaround = 1;

% time stretching settings
Stretch.step = 0.01;
Stretch.nbef = 1.5/Stretch.step;
Stretch.naft = 1.5/Stretch.step;
Stretch.gkern = 0.05;

% Realtime settings
Realtime.gkern = 0.05;
Realtime.step = 0.01;
Realtime.nbef = 1.5/Stretch.step;
Realtime.naft = 1.5/Stretch.step;
Realtime.nbet = 0.5/Realtime.step;
Realtime.naround = 1/Realtime.step;

% number of bootstrap sims to run
nboot = 10000;

Regression.gkern = 0.015;
% Time to look beyond starts and stops for regression periods
Regression.reachpad = 0.25;
% Time to look into hold period for error calculations
Regression.errorwindow = 1;
Regression.ratesuse = 'slidenormrates';
% Regression parts, no error
Regression.parts = {'Pos', 'Vel', 'Dir', 'Speed', 'Acc', 'AMag', 'ASign'};
Regression.parts_long = {'3D Position', '3D Velocity', '3D Direction', 'Speed', '3D Acceleration', 'Accel. Magnitude', 'Signed Accel.'};
% Regression parts, error
Regression.parts_err = {'Pos', 'Vel', 'Dir', 'Speed', 'Acc', 'AMag', 'ASign', 'Err', 'EMag'};
Regression.parts_err_long = {'3D Position', '3D Velocity', '3D Direction', 'Speed', '3D Acceleration', 'Accel. Magnitude', 'Signed Accel.', '3D Error', 'Error Magnitude'};
% Set lag parameters
Regression.maxlagsecs = 1;
Regression.maxlagsamps = Regression.maxlagsecs*vidframerate;
% segmentation types for regression
Regression.segstr = {'Full', 'Reach', 'ReachE', 'Hold'};
Regression.cutoff = 6;

% which smooth for plotting
Plot.gkern = 0.05;

save([savedir 'SS.mat']);


% MISSING FILE NOTES
% VM01 and VM04 were dystonic tremor and are bad / weird. Don't use
% VM11 has no units, also only 1 session
% VM13 was aborted due to shoulder pain
% VM15 behavior was essentially unusable
% VM24 didn't happen
% VM25_2 has bad anipose