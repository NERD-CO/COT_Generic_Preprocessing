Rex Tien 1/21/24

Preprocessing code for center-out task DBS data [NERD CO]

This readme walks through the steps of preprocessing for center-out task DBS data. There are several steps, outlined below, and some steps require processing with separate programs outside MATLAB.

The code relies on several environmental variables to set your working directories. Set these and then restart MATLAB before running.
'SYNOLOGYDIR': The path to the "CenterOutTask" folder in the rsrch5 synology drive ( something like '\\som-nsg-rsrch5.ucdenver.pvt\documents\CenterOutTask')
'DBSDATADIR': The path to your general dbs data folder. Each subject will get a subfolder in this directory with the raw data files.
'AOCONVERTERDIR': The path to your install of the AlphaOmega mpx file converter software (something like 'C:\Program Files (x86)\AlphaOmega\Converter')
'ANACONDAPYDIR': The path to your anaconda python install (something like 'C:\Users\username\anaconda3\envs\DEEPLABCUT\python.exe')
'PLXSDKDIR': (OPTIONAL) The path to where you downloaded the Plexon MATLAB offline SDK (optional, if you want to do hand-sorting in plexon)
'PLXDIR': (OPTIONAL) The path to where you want to save the plx files, where you will do sorting
'DLCDIR': The path to where your DLC projects are stored

It also requires the following outside packages:
AlphaOmega MPX file converter. Windows x64 installer is included in this repo
DeepLabCut, installed with Anaconda. Installation instructions here: https://deeplabcut.github.io/DeepLabCut/docs/installation.html
Plexon SDK (if you want to do Plexon OfflineSorter sorting. Can be found here: https://plexon.com/software-downloads/#software-downloads-SDKs). Download the MATLAB Offline Files SDK.zip, add it to your MATLAB path)

And requires the following files in your 'DBSDATADIR' directory (files have been provided here, move them to DBSDATADIR:
COT_CamAlign_Fixes.xlsx
COT_ReachFind_Fixes.xlsx

1. Run PullAndPrep: this creates the subject folder in 'DBSDATADIR' and does initial Preprocessing. It runs the following subroutines:
	a. PullFromSynology - pull the video, AO data, task data and notes file from Synology, convert AO mpx files to mat files
	b. AutoGenKey - scans the available files and tries to automatically generate a key file saying which AO files go with which Task and Video files
	c. Batch_AO2plx - (optional) batch-converts AO files to plx files for hand-sorting with Plexon OfflineSorter
	d. GenRunDLC - creates a DLC project folder with good and runs the k-means clustering for labeling
	
2. Now run the steps to complete the DLC processing in the new folder that was created in your 'DLCDIR'\SS##_ndimD directory
	- Run SS##_label.py, label the frames
	- Run SS##_train.py, wait
	- Run SS##_eval-analyze.py, wait

3. (OPTIONAL) Sort units in 'PLXDIR'\SS## folder, save sorted files in 'PLXDIR'\SS##\Sorted

4. Run CheckAlignment to examine TTL and camera frame alignment. If alignment isn't right, make changes to COT_CamAlign_Fixes.xlsx to fix the issue

5. Run AnnotateVids to create annotated videos in your 'DBSDATADIR'\SS##\Annotated_Video\ folder
	- Examine these videos, edit the accompanying .xlsx annotations files to exclude bad trials
	- ******** You must enter a 'g' in cell D1 to indicate that you have reviewed the file
	
6. (OPTIONAL) creates and run a SetSS version to generate your settings for subject codes, preprocessing in the next step

7. Run Preprocess_All: this does the steps to generate a file with reach events, kinematics and neural data. Note that it should be run first with 'doplot' set to true, to inspect reach event finding. It runs the subroutines:
	a. extract_dlc_time_and_pos
	b. fix_frameskips
	c. fix_blips_and_gaps
	d. calculate_kinematics
	e. calculate_ktime_rates
	f. detect_reach_events
****** It is highly recommended to run with 'doplot' set to true first, and to step through and examine all reach event detection plots. If strange detection is identified, either screen the trial by putting an 'x' in the video annotation .xlsx, or edit the COT_ReachFind_Fixes.xlsx file to make detection better.

Preprocess_All will output Reach_Kin_Neu_SS##_##.mat files with the following fields:
FR: Firing rates (FRs) of the units that that were hand sorted
	SmoothG15: FRs smoothed with a 15ms Gaussian
	SmoothG50: FRs smoothed with a 50ms Gaussian
Kin: Kinematics of fingertip (at camera frame times, note that NaNs may be present due to tracking / frame skip gaps)
	Raw: Kinematics without smoothing
		Pos: position
		Vel: velocity
		Speed: magnitude of velocity
		Dirs: direction of velocity
		Acc: acceleration
		AMag: magnitude of acceleration
		ASign: component of acceleration in the same direction as velocity
	SmoothG15: same as "Raw" but smoothed with a 15ms Gaussian
	SmoothG50: same as "Raw" but smoothed with a 50ms Gaussian
N: Neural data
	Depth: Depth of recording above planned target
	SpkID: channel info for hand-sorted units (note: does not exist if you didn't do any hand-sorting)
	SpkTimes: spike times for hand-sorted units (note: does not exist if you didn't do any hand-sorting)
	SPK: CSPK_## fields from the AlphaOmega data
	LFP: CLFP_## fields from the AlphaOmega data
Reach: Reach event timestamps / reach info. Note: fields entries only exist for "valid" reaches
	reachdirs: directions of reach (degrees, 0 degrees is right, counterclockwise)
	reachnum: original trial number of the reach
	outreach: true if reach was center-out, false if reach was out-center
	reachpeaks_rough: initial peak speed in each reach used for reach timing finding
	reachpeaks: time of peak speed in each reach, in the reach direction
	reachstarts: time of each reach start
	reachstops: time of each reach stop
	reactiontime: time between go cue and reach start
Time: Time info
	ktime: timestamps for the kinematics fields
	ArrowIsShown: times when the arrow cue is shown
	CenterIsShown: times when the center target is shown (go cue for out-center reaches)
	StartOfMovement: times when the center target is removed and outer target is shown (go cue for center-out reaches)