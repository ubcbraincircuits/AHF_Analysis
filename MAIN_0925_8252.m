%  Script FARM_toTim_reps
%  Script that calls functions together in order to generate the following:
%     FARM_Preprocess : Function who load the data on the workspace
%  *  FARM_Process : Function who process the data before exploration
%     FARM_seed : Function who seed the calcium activity and compare it to behavior
%     FARM_video : Function who play or record video of dual calcium and behavior
%     FARM_averaging : Function who average trial and seed

%% Open a pair of calcium and behavior data
clear all
close all, clc

% varOpen would be the structure array that acts as the working folder where most data is kept, such as the parameters that are set for each of
% the functions that are being run.  

% CAGE GENERATION 1 ----------------------------------------------------------
%varOpen.file = 'D:\Ai94_10_mouse_cage_march_2018_males\20180518\3\Videos\'; 
%varOpen.filesave = ['D:\Ai94_10_mouse_cage_march_2018_males\20180518\test']; 
%varOpen.filecam =  'D:\Ai94_10_mouse_cage_march_2018_males\videos_eye_may_2018\';
%varOpen.filecam2 = 'D:\Ai94_10_mouse_cage_march_2018_males\videos_bottom_april_2018_auto\';
%varOpen.file = 'E:\07302018_Tet_\20180825\3\Videos\'; 

varOpen.working_folder = 'D:\20180925-252\Processed\';

varOpen.file(1).rep = 'D:\07302018_Tet_\20180925\3\Videos\'; % load the text files
varOpen.file(2).rep = 'D:\07302018_Tet_\20180926\3\Videos\'; 
varOpen.file(3).rep = 'D:\07302018_Tet_\20180927\3\Videos\';
varOpen.file(4).rep = 'D:\07302018_Tet_\20180928\3\Videos\';
varOpen.file(5).rep = 'D:\07302018_Tet_\20180929\3\Videos\';
% varOpen.file(6).rep = 'D:\07302018_Tet_\20180929\3\Videos\';
%varOpen.file(4).rep = 'E:\07302018_Tet_\20180917\3\Videos\';
%varOpen.file(5).rep = 'E:\07302018_Tet_\20180828\3\Videos\';
%varOpen.file(6).rep = 'E:\07302018_Tet_\20180829\3\Videos\';
%varOpen.file = []

%%varOpen.file(1).rep = 'I:\Cage1\20180225\3\Videos\'; varOpen.file(2).rep = 'I:\Cage1\20180225\3\Videos\'; 
%%varOpen.file = [];

% for loading preprocessed files the beginning of filename must shown in filesave without mouse #
varOpen.filesave = 'D:\20180925-252\Processed\20180925_29'; % load the videos 
varOpen.filecam =  'D:\07302018_Tet_\video_eye\'; % where behavioural videos are stored (raw data) 
varOpen.filecam2 = 'D:\07302018_Tet_\video_front\';
varOpen.filecam3 = 'D:\07302018_Tet_\video_bottom\';
%
varOpen.mouse = '8252'; 
% varOpen.mouse = '8474'; 
% varOpen.mouse = '8423'; % name of mouse
% Use 0 or 1, 0 will load previously processed data from .mat files, 1 will allow re-processing of previously processed raw files.
varOpen.reprocess = 0; 
%file = 'D:\Ai94_10_mouse_cage_march_2018_males\20180518\3\Videos\';
%filecam = 'D:\Ai94_10_mouse_cage_march_2018_males\videos_eye_may_2018\';
%filecamSide = 'D:\Ai94_10_mouse_cage_march_2018_males\videos_bottom_april_2018_auto\';
%filesave = 'D:\Ai94_10_mouse_cage_march_2018_males\20180518\'; % saving dataset in these files (to reopen faster

% CAGE GENERATION 2 ----------------------------------------------------------
%%varOpen.file = 'I:\Cage1\20180225\3\Videos\M2016090793_stim_1519623845.raw';
%%varOpen.file(1).rep = 'I:\Cage1\20180225\3\Videos\M2016090793_stim_1519623845.raw'; varOpen.file(2).rep = 'I:\Cage1\20180225\3\Videos\M2016090793_stim_1519623845.raw';
%%varOpen.file(1).rep = 'I:\Cage1\20180225\3\Videos\'; varOpen.file(2).rep = 'I:\Cage1\20180225\3\Videos\'; 
%%varOpen.file = [];
% varOpen.file = 'I:\Cage1\20180225\3\Videos\'; 
% varOpen.filesave = ['I:\AHF\dataset_Cage1_180225_']; 
% varOpen.filecam =  'I:\Cage1\front\';
% varOpen.filecam2 = 'I:\Cage1\side\';
% varOpen.mouse = '793'; 
% varOpen.reprocess = 1;

% STROBING EXAMPLE ----------------------------------------------------------
% varOpen.file = 'D:\CBV\2017-02-22\2017_02_22_Ai94_2.tif'; 
% varOpen.filecam = 'D:\CBV\2017-02-22\2017_02_22_Ai94_2.h264'; 
% varOpen.SF = 150;
% varOpen.reprocess = 1;
% varOpen.filesave = 'D:\CBV\2017-02-22\2017_02_22_Ai94_2.mat'; 

% YUKI'S EXAMPLE ----------------------------------------------------------
% varOpen.file = 'E:\TeamCa2\Y1\Y1_Nov21pm_15Hz_8x8.tif';
% varOpen.filecam = 'E:\TeamCa2\Y1\Y1pm_2014-11-21 10-53-21.810.wmv';
% varOpen.SF = 15;
% varOpen.filesave = 'E:\TeamCa2\Y1\test'; 
varOpen.L=8.2; % width of image: tif= 8.2, raw= 10.25 
varOpen.delete = 0 % whether you want to delete the temporary behavioural and averaging files 1= delete 0= no 
[expe,varOpen] = FARM_PreProcess_teri(varOpen); % expe contain data like the date, file name, bregma pixels, etc. of each of the images

%% preprocessing
%expe([39])=[]
% expe([151])=[] %exclude trial 151 
% remove images with unwanted motion by looking at the video displayed in this section. 
varOpen.varProc.run_list=[1:4, 9:25, 27:33, 35:38, 40:41, 43, 45, 49:54]; % Indicate the expe files you want displayed. 
%expe([59])=[]
%varOpen.varProc.run_list=[57];
varOpen.varProc.tb = 3; % temporal binning- 1 indicates no change
varOpen.varProc.useBe1=0; % different behavioural cams to monitor different parts of the body 
varOpen.varProc.useBe2=1; % 1 = turn on; 0= turn off 
varOpen.varProc.useBe3=1;
varOpen.varProc.Ca_flipLR=1; % flip image of calcium data(sometimes displayed other way around)
%varOpen.varProc.L=8.2;
varOpen.blueCorr=1; % reflectance correction: 0= no, 1= yes
%varOpen.varProc.delete=1; % delete temp files for previous behavioural and averaging data; 1= yes, 0= no ((double check with preprocess teri)) 
[expe,varOpen] = FARM_Process_beh(expe,varOpen); 

% Change the parameters shown above according to changes that need to be made(shifts and glitches in the video, reflecting the image, behaviour
% cams) after viewing the video

% The following code is used to write an avi file of the Reflectance, Calcium and Behaviour videos.  

% Check to see if there are more than 1 trial
if size(varOpen.varOutput.REFv,3)>1 
    % Will display a video if there is more than one trial  
    % OIA_new_video function accepts 2nd arg called v which has options for saving
    % v.video file is the filename to save
    v.videofile=[varOpen.working_folder varOpen.mouse ];
    
% v.mv: set this to the value corresponding to what you would like to view 
%           =1 if stored in a avi
%           =2 if stored in a mat 
%           =3 if just play the video in loop
%           =4 if just the output

% If options 1 or 2 were chosen, each of the trials will be displayed and will advance every time you click a key.  Use this to check which trials
% should be used for further analysis. 

% If the video is saved as an avi, the avi file can be opened in ImageJ as well and further accessed to see which images are acceptable for
% further analysis.  Adjust the parameters above in the run list and rerun this section if images are to be excluded. The video will only be played
% once. 

    v.mv=1;
    [~,CC]=OIA_new_video(varOpen.varOutput.REFv,v);

    figure; % This creates the video combining all the frames together 
    for i=1:size(CC,4)
        imagesc(CC(:,:,:,i))
        title(num2str(i))
        pause
    end
    
% If it only contains one trial, then it will display the image 

else
    imshow(varOpen.varOutput.REFv); colormap gray;
end

%% Seed region body vs calcium
clear CC i v
varOpen.Seed.n=3 % number of rois to select
% if this section has already been run, it will not prompt you to draw the roi again. If you wish to observe another ROI, choose the redraw option.
varOpen.Seed.redo_roi=0; % roi need to be redrawn? 0= no, 1= yes
% (ex. leaving it [] will show you the data for the whole trial) 
varOpen.Seed.win = [] % window being displayed in seconds([min, max]) 
% For the colour bar option, the command window will display the min and max values for correlation. If you want to adjust the scaling yourself,
% enter the values for scaling you want in the form [min max].  If not,putting [] will have the function set its own scaling. 
% The command window will display the minimum and maximum values to give you an idea of the scaling that should be used.  You can adjust colour
% bar according to those values. 
varOpen.Seed.colour_bar = [] % indicates the scaling for the colour bar (third image on the first row) 
[varOpen,expe] = working_seed_teri(varOpen,expe);

% RIGHT CLICK WHEN DONE SELECTING THE CORRECT PIXEL, DOING THIS WILL SAVE
% BEHAVIOUR DATA PROPERLY FOR AVERAGING! 

% Choose the seed you would like to look into and see the correlation between the movements and the seed chosen.
% The image on the very left indicates the correlation between the calcium activity of the seed and the real-time movement of the mouse associated with it 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYBRID CORRELATION (FIRST IMAGE ON THE FIRST ROW) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shows the correlation between the calcium data of the seed chosen and the movement of the mouse that correlates to it. 
% (i.e closer it is to red, the higher connectivity between the movement in the area and the seed chosen)  

%%%%%%%%%%%%%%%%%%%%% REFLECTANCE AND SEED CORRELATION (SECOND IMAGE ON THE FIRST ROW) %%%%%%%%%%%%%%%%%%%%%
% Reflectance image is just there to help with visualization
% The seed image shows the correlation between the seed chosen and the other pixels in the brain 
% The white outline indicates the division of areas in the brain according to the Brain Allen Reference Atlas 

%%%%%%%%%%%%%%%%%%%%%%%%%% HYBRID CROSS-CORRELATION (THIRD IMAGE ON THE FIRST ROW) %%%%%%%%%%%%%%%%%%%%%%%%%%
% Shows calcium activity in the areas of the brain which are responsible for movement in the ROI chosen.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE PLOT ON THE BOTTOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cyan dots: instantaneous rates of licks
% magenta line: curve fitting line between the instantaneous rates of licks.  
% red dots: time in which licks occurred 
% blue asterisk: indicates when a reward was given 
% green triangles: cue for the start of a go task event (lick)
% red triangles: cue of the start of a no go task event (no lick)
% varying coloured lines in the middle: (behavioural) movement of the different ROIs
    % movement tracked via 3 behavioural cameras and are all syncronized;
    % can compare it with the fluorescence data located below it
% Refer to legend to identify fluo, refl and sub data. This data is gathered via brain imaging. 
% Red vertical lines are evenly spaced out to indicate when the trial ends and the next one starts(data is concatenated) 

%% visu video and/  profile at the same time

varOpen.Video.positionWin = [4750 4783]; % position (in seconds)
% The range indicated below shows the DF/F that should be displayed in the animation. This is NOT regarding the range of DF/F you want to see in the plot on the right.  
varOpen.Video.range = [-10 10]; % range of DF/F (%) displayed  
varOpen.Video.description = 1 % show (1) or don't show (0) descriptions of the trials at the top of the plot (like the period, wait time, etc) 
varOpen.Video.rec = 0; % record AVI? 0=no (just watch), 1=yes (only video), 2= all (with the graph)
varOpen.Video.videofile = 'E:\07302018_Tet_\20180925\3\Videos\0925_29_8423_fluo_refl.avi'; % videofile name

FARM_video_teri(varOpen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE PLOT ON THE LEFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Black dot would mark the point you selected in the previous section (FARM_seed). 
% The animation on the left would be the animation that corresponds to the fluorescence. With the proper range set, you will be able to visualize the change in the animation
% in accordance to the data plotted in the graph on the right.  
% When overall fluorescence is high, the image approaches red.  When fluorence is low, it approaches blue. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE PLOT ON THE RIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% red dots: time stamps when the mouse licked. 
% blue asterix: shows a reward received
% green triangle: cue for the start of a go task event (lick)
%   GO TRIALS:
%         +2: reward
%         -2: missed 
%         -4: licked early 
%   NO GO TRIALS
%         +1 withheld licking
%         -1 licked in the correct interval but wrong action
%         -3 licked early, no lick 
% red triangle: cue for the start of a no go task event (no lick) 


% The fluorescence and reflectance is shown on the bottom. The subtraction line is shown as well, but it is shifted upwards in order to 
% clearly visualize its corrected calcium activity in relation to its movement (added difference to the max value of calcium data shown in FARM_video)

% Colours of marked ROIs correspond to the colour of the lines displayed in the plot. 
% (The L- looking line) The black line horizontal line indicates a second of the trial and the black vertical line represents 10% of the fluorescent data. 

%% event averaging

% Will display a montage of brain images as well as its corresponding time,
% the average brain image of max and min activity,
% plot concerning fluo and refl data and the behavioural gradient of all ROIs 

varOpen.ave.mvt_roi = 3; % The ROI that should be used to take kSD into account regarding movement threshold
varOpen.ave.mvt_kSD = 4; % will discard trial if mvt above threshold (mvt_kSD x std dev)
varOpen.ave.mvt_pre = 0.5; % check mvt for which regions from range 0 to 1 0.5: check first half of record; 1: check the whole window 

varOpen.ave.ntrials = [] % select the trials you want to analyze (leave as [] if you want to see all trials) 
varOpen.ave.redo = 1 % restart the whole process (change if any parameter except for n trials is changed) ? 0=no, 1=yes
% q-type: -1= behavioural ROI , 0:lick, 1:reward given , 2:GO=2, 3:GO=-4, 4:GO=1, 5:GO=-1, 6:GO=-3; 7:GO=-2 GO refers to its response to stimulus 
qtype_list= [3]; % specify qtype corresponding to GO trial stated above 
varOpen.ave.pre = 2.5; % Time in seconds analyzed before and after the cue 
varOpen.ave.refr_period = 3; % refractory period: period with no response after stimulation
varOpen.ave.temp_filter = 0.1; % temporal filter: smooths out the frequency 
varOpen.ave.range = [-4 10]; % Range of DF/F (change in fluorescence) in % 
varOpen.ave.subtraction_line = 1; % If a subtraction(for fluorescence and reflectance) line is wanted, set subtraction_line = 1, if not, set it to 0.
varOpen.ave.seed_pixel = [52,27]; % Pixel dim. of seed you would like to analyze ([x,y]).  Leave it as [] for the interactive plot to click on various regions of the brain to analyze its activity

qtype_files= {};
for i= 1:length(qtype_list) 
    qtype_files{i} = ['qtype' num2str(qtype_list(i)) '.mat'],
end 

varOpen.ave.qtype_folder = varOpen.working_folder;
for m= 1:length(qtype_list)
    try
    varOpen.ave.qtype = qtype_list(m);
    varOpen.ave.outfile=[varOpen.ave.qtype_folder qtype_files{m}];
    [varOpen,expe] = FARM_averaging_FINAL(varOpen,expe);
    catch
        disp('Done with interactive plot, continuing loop');
    end
end

clear i m 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% red dotted vertical lines: time corresponding to the brain montages above.  It shows you the activity of the brain at that specified time. 
% green line: fluorescence 
% blue line: represent the reflectance
% black line: subtraction of reflectance from fluorescence
% different colours of ROI: related to movement of each ROI 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 2 (Will only pop up when you decide to redo the process) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% blue triangles: depends on the q-type chosen 
    % for q-types:
        % 0: the time when the mouse licked 
        % 1: time when reward was given 
        % others: time the GO trial began 
% black line: represents the filtered behaviour of the mouse (filtered behaviour)
% red line: the spikes show the time stamps in which the calcium levels exceeded the threshold 

%% Single Trial Data- Loading it 
% This function will create brain montages of brain activity with the progression of time. The frames that are being taken into account are
% specified in the parameters below. The graphs beside the brain montages correspond to the movement of the mouse.  Every new row represents a new
% trial. 

% Figure one would show all the brain montages and the movements of the mouse in the specified trials.  The second figure would take the
% average of the brain images and movement from all trials, generating average brain montages and movement.

% These files are automatically saved to the working folder specified in preprocessing (folder names are the qtypes of the trials). 
varOpen.singletrial.trial_interval = [400:423] % state the interval of trials (50 trials would be the max no.) you want to view ([min:max]) 
varOpen.singletrial.time = [1:50] % frame interval being observed [min:max] or [] for whole time interval
varOpen.singletrial.temp_bin= 50 % n of frames you want to display in the montage; make sure it is a multiple of the # of seconds shown in singletrial.time
varOpen.singletrial.qtype_list = qtype_list; % making qtype_list a field so it can be a part of varOpen and be executed in FARM_singletrial 
varOpen.singletrial.qtypeoutfile = varOpen.working_folder; 

FARM_singletrial(varOpen) 
