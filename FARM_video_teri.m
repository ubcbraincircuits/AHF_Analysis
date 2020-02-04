% function FARM_video(varOpen)
% Function who play or record video of dual calcium and behavior
%    FARM_Preprocess : Function who load the data on the workspace
%    FARM_Process : Function who process the data before exploration
%    FARM_seed : Function who seed the calcium activity and compare it to behavior
%  * FARM_video : Function who play or record video of dual calcium and behavior
%    FARM_averaging : Function who average trial and seed
%
% INPUT: varOpen.Video. (Structure)
%    .rec: (scalar) recording mode (0: just display [default], 1: records just image, 2: records image and graph)
%    .fps: (scalar, in Hz) frame per second [default: 10]
%    .positionWin: (matrix of 2 scalars, in sec) time position [start end] or [] for all (default)
%    .range: (matrix of 2 scalars, in DF/F%) range displayed [min max] or [] for full range (default) 
%    .videofile: (string) avi file name (default: [local_folder '\temp.avi'])
%
%         NB: to remove any field of a structure, use rmfield:
%             e.g. varOpen = rmfield(varOpen,'Video'); to use only default
%             e.g. varOpenv.Video= rmfield(varOpenv.Video,'fps'); to only remove one sub-field

function FARM_video(varOpen)

% --------------  1. define default variable ------------------------------
disp('---------------------------------------------------------------------')
try, disp(['              varOpen.Video.rec=' num2str(varOpen.Video.rec)]); 
    catch, varOpen.Video.rec = 0; disp('(not defined) varOpen.Video.rec=0'); end; % set the default value of video display = 0 (=no recording)   
try, disp(['              varOpen.Video.fps=' num2str(varOpen.Video.fps)]); 
    catch, varOpen.Video.fps = 10; disp('(not defined) varOpen.Video.fps=10'); end % set the default display frequency at 10Hz
try, disp(['              varOpen.Video.positionWin=' num2str(varOpen.Video.positionWin)]); 
    catch, varOpen.Video.positionWin = [0 size(varOpen.varOutput.I,3)/varOpen.varOutput.SFr]; % set the default time window at full sequence (all frames)
        disp(['(not defined) varOpen.Video.positionWin = ' num2str(varOpen.Video.positionWin)]); end
if length(varOpen.Video.positionWin)==0, varOpen.Video.positionWin=[0 round(size(varOpen.varOutput.I,3))/varOpen.varOutput.SFr]; end     
try, disp(['              varOpen.Video.range=' num2str(varOpen.Video.range)]); % set the default range display to min and max of the data
    catch, varOpen.Video.range = [min(varOpen.varOutput.I(:)) max(varOpen.varOutput.I(:))]; 
        disp(['(not defined) varOpen.Video.range = ' num2str(varOpen.Video.range)]); end
try, disp(['              varOpen.Video.videofile=' varOpen.Video.videofile]); % set the default videofile name in the local directerory as temp.avi
    catch, varOpen.Video.videofile = [pwd '\temp.avi']; 
        disp('(not defined) varOpen.Video.videofile= temp.avi file in local folder'); end   
try, disp(['              varOpen.Seed.blue_correction=' num2str(varOpen.Seed.blue_correction)]); % set the hemodynamic correction on (=1)
    catch, varOpen.Seed.blue_correction = 1; disp('(not defined) varOpen.Seed.blue_correction=1'); end 
[task_sound,Fs] = audioread('task.wav'); % load the sound "task" and "drop"
[drop_sound,Fs] = audioread('drop.wav');
position = round(varOpen.Video.positionWin*(varOpen.varOutput.SFr))+1; % convert time window (in sec) in frames
duration_video = (position(2)-position(1))/varOpen.Video.fps; % calculate the total duraction of the video
sound_band = zeros(round(duration_video*Fs),1,'single'); % generate the vector of sound  (wav file)
lick_sound = sin(880*2*pi*[1:round(Fs/20)]'/Fs); % generate the sound for "lick" (= a simple A beep at 880Hz)  
delta_sound=length(sound_band)./(position(2)-position(1)); % ratio of sound tape relatively to data

seed_temp_file = [varOpen.folder 'seed_temp.mat']; % file where the last seeding data were stored for display (See FARM_seed.m)
%seed_temp_file = [varOpen.folder 'Ave_temp.mat']; % file where the last seeding data were stored for display (See FARM_seed.m)
load(seed_temp_file);
try, seed_info.re; catch, varOpen.Seed.blue_correction = 0; end % if no reflectance channel present, set the reference correction to zero



% --------------  2. display the background _------------------------------
fig1 = figure;
fig1.Renderer = 'Painters';
subplot(122)
plot(seed_info.bxx,seed_info.ca,'g','LineWidth',2); hold on; % display calcium and eventually reflectance
try, plot(seed_info.bxx,seed_info.re,'b','LineWidth',2); end

try, plot(seed_info.bxx,max(seed_info.ca)+(seed_info.ca-seed_info.re),'k','LineWidth',2); end % display eventually calcium corrected (substracted) by reflectance
plot([varOpen.Video.positionWin(1)+.5 varOpen.Video.positionWin(1)+1.5],[min(seed_info.ca) min(seed_info.ca)],'k','LineWidth',2); % display the x-scale 1 sec
 text(double(varOpen.Video.positionWin(1))+1.5,double(min(seed_info.ca)),'1sec')
plot([varOpen.Video.positionWin(1)+.5 varOpen.Video.positionWin(1)+.5],[min(seed_info.ca) min(seed_info.ca)+10],'k','LineWidth',2); % display the y-scale 10%
 text(double(varOpen.Video.positionWin(1)),double(min(seed_info.ca)+10),'10%')

if varOpen.varOutput.sBe ~= 0 % if there is behavior data
    ben = OIA_n(seed_info.be) * (max(seed_info.ca)-min(seed_info.ca)); % normalize the behavior movement in the same range of calcium to make comparison easier
    for i = 1:size(ben,1)
        plot(seed_info.bxx,max(seed_info.ca)+ben(i,:),'Color',seed_info.couleur(i,:),'LineWidth',2); % display each behavior from each body roi
    end
end
plot(varOpen.varOutput.event_list,seed_info.base_beh.*ones(size(varOpen.varOutput.event_list)),'r.'); % eventually, display each lick event
plot(varOpen.varOutput.event_list_reward,seed_info.base_beh.*ones(size(varOpen.varOutput.event_list_reward)),'b*'); % eventually, display each task event

for i= 1:length(varOpen.varOutput.event_list_task_info)
        if (length(strfind(char(varOpen.varOutput.event_list_task_info(i)),'GO=2'))|length(strfind(char(varOpen.varOutput.event_list_task_info(i)),'GO=-2'))|length(strfind(char(varOpen.varOutput.event_list_task_info(i)),'GO=-4'))) ==1
            plot(varOpen.varOutput.event_list_task(i),seed_info.base_beh,'gv');
        else 
            plot(varOpen.varOutput.event_list_task(i),seed_info.base_beh, 'rv'); 
        end
end

% plot(varOpen.varOutput.event_list_task,seed_info.base_beh.*ones(size(varOpen.varOutput.event_list_task)),'gv'); % eventually, display each drop delivered

for i = 1:length(varOpen.varOutput.event_list_task_info) % for each task event, write the text information of the task
    if ((varOpen.varOutput.event_list_task(i)>varOpen.Video.positionWin(1))&(varOpen.varOutput.event_list_task(i)<varOpen.Video.positionWin(2)))
        str_event_task = char(varOpen.varOutput.event_list_task_info(i));
        str_event_task_list_1= str_event_task(1:30);
        str_event_task_list_2= str_event_task(31: length(str_event_task));
        if varOpen.Video.description == 1 
         text(double(varOpen.varOutput.event_list_task(i)),double(seed_info.base_beh)-5,str_event_task_list_1,'Rotation',90,'FontSize',8);
         text(double(varOpen.varOutput.event_list_task(i))+0.5,double(seed_info.base_beh)-5,str_event_task_list_2,'Rotation',90,'FontSize',8);
        else 
        end 
    end
end

for i = 1:length(varOpen.varOutput.jonction) % for each jonction, the position of each concatenation between different trials
    position_jonction = sum(varOpen.varOutput.jonction(1:i))/varOpen.varOutput.SFr; % convert the jonction location in sec in real frame
    plot([position_jonction position_jonction],[min(seed_info.ca) seed_info.base_beh],'r')
end

xlim(varOpen.Video.positionWin) % restrict the display to the window
xlabel('time (s)'); ylabel('DF/F (%)') % display the axis
if varOpen.Seed.blue_correction == 1, legend('fluo','refl','subs', 'location', 'southeast'); end % legend only if reflectance present
hold off



% --------------- 3. prepare the vector of events -------------------------
JJL = zeros(size(seed_info.ca)); % prepare the vector of events: lick, task and drop
JJT = zeros(size(seed_info.ca));
JJR = zeros(size(seed_info.ca));
for i = 1:length(varOpen.varOutput.event_list) % for each lick
    pos = 1+round(varOpen.varOutput.event_list(i).*varOpen.varOutput.SFr); % establish the location (in frame)
    JJL(pos) = 1; % add an event at the calcium location
end;
for i = 1:length(varOpen.varOutput.event_list_reward)
    pos = 1+round(varOpen.varOutput.event_list_reward(i).*varOpen.varOutput.SFr);
    JJR(pos) = 1;
end;
for i = 1:length(varOpen.varOutput.event_list_task)
    pos = 1+round(varOpen.varOutput.event_list_task(i).*varOpen.varOutput.SFr);
    JJT(pos) = 1;
end;
JJL = uint8(JJL); % convert in 8bits (for video-image)
JJR = uint8(JJR);
JJT = uint8(JJT);



% --------------- 4. loop of display/recording of the data ----------------

if varOpen.Video.rec ~= 0, % if video recording
    max_ite = 1; % set only one loop
    vidObj = VideoWriter(varOpen.Video.videofile); % open the file to write
    vidObj.FrameRate = varOpen.Video.fps; % set the framerate
    open(vidObj);
else
    max_ite = 1000000; % set a large number of loop
end

LUT = imresize([100:-1:0]',[size(varOpen.varOutput.I,1)/2 size(varOpen.varOutput.I,2)/10 ]); % generate the LUT on the edge
LUT = OIA_color(LUT,jet);

for ite = 1:max_ite % for each iteration
for i = position(1):position(2) % for each location (in frame)
    
    
    
    % --------- 4.1 generate the images IMG -------------------------------
    IMG1 = varOpen.varOutput.I(:,:,i,1); % pickup the calcium data of this specific frame
    try, IMG2 = varOpen.varOutput.I(:,:,i,2); end % and eventually, the reflectance frame
    if varOpen.Seed.blue_correction == 1, IMG = IMG1 - IMG2; else, IMG = IMG1; end % transfert substaction of fluo only to IMG (if not reflectance)

    if varOpen.varOutput.sBe ~= 0 % if there is mvt data
        BEH = single(varOpen.varOutput.Be(:,:,i)); % pickup data of behavior
        IMG = imresize(IMG,size(BEH,2)/size(IMG,2),'nearest'); % scale IMG at the same size than BEH
        BEH = OIA_color(BEH,gray); % convert in RGB
        for ii = 1:size(seed_info.roiBeh,3) % for each behavior roi
            roiBeh_disp = seed_info.roiBeh(:,:,ii); % pick up the roi
            roiBeh_disp = imdilate(roiBeh_disp,ones(3)) - roiBeh_disp; % generate the contour and imprint if the the behavior video
            for j = 1:3, 
                BEH(:,:,j) = BEH(:,:,j).*uint8(~roiBeh_disp) + uint8(255*roiBeh_disp.*seed_info.couleur(ii,j)); 
            end;
        end
    else
        IMG = imresize(IMG,3,'nearest'); % no behavior, expand the calcium data so it is easier to see
    end
    
    % cutoff the data to the range:
    IMG = IMG.*((IMG>varOpen.Video.range(1))&(IMG<varOpen.Video.range(2))) + varOpen.Video.range(2).*(IMG>=varOpen.Video.range(2)) + varOpen.Video.range(1).*(IMG<=varOpen.Video.range(1));
    IMG(1,1) = varOpen.Video.range(1); IMG(1,2) = varOpen.Video.range(2);
    IMG = OIA_color(IMG,jet); % convert in RGB
    ratio = size(IMG,1)/size(varOpen.varOutput.I,1); 
    x3d = round(seed_info.x3*ratio); % inprint the location of the reference seed to the image (with a cross) 
    y3d = round(seed_info.y3*ratio);
    IMG(y3d-10:y3d+10,x3d,:)=zeros(21,1,3,'uint8'); 
    IMG(y3d,x3d-10:x3d+10,:)=zeros(1,21,3,'uint8');
    IMG(y3d-4:y3d+4,x3d-4:x3d+4,:)=255*ones(9,9,3,'uint8'); 
    IMG(y3d-3:y3d+3,x3d-3:x3d+3,:)=zeros(7,7,3,'uint8');    
    for j = 1:3, IMG(:,:,j) = IMG(:,:,j).*uint8(imresize(varOpen.varOutput.roi,[size(IMG,1) size(IMG,2)],'nearest')); end; % add the roi mask
    
    IMG(1:size(LUT,1),1:size(LUT,2),:)=LUT; % add the LUT on the top left
     IMG = insertText(IMG,[0 0],num2str(varOpen.Video.range(2)),'Font', 'Arial','FontSize',10,...
                         'BoxColor', 'red', 'BoxOpacity', 0,'TextColor','white'); % write the maximum at the top of the LUT
     IMG = insertText(IMG,[0 size(LUT,1)*.8],num2str(varOpen.Video.range(1)),'Font', 'Arial','FontSize',10,...
                         'BoxColor', 'red','BoxOpacity',0,'TextColor','white'); % write the minimum at the top of the LUT
  
    BEHx = round(size(IMG,2)/20); 
    BEH(1:BEHx,1:BEHx,1)=255*JJL(i).*ones(BEHx,BEHx,1,'uint8'); % if JJL (See before) = 1, put a red square in the top left
    BEH(1:BEHx,1:BEHx,2:3)=zeros(BEHx,BEHx,2,'uint8');

    BEH(1:BEHx,BEHx+[1:BEHx],3)=255*JJR(i).*ones(BEHx,BEHx,1,'uint8'); % same for JJR
    BEH(1:BEHx,BEHx+[1:BEHx],1:2)=zeros(BEHx,BEHx,2,'uint8');

    BEH(1:BEHx,2*BEHx+[1:BEHx],2)=255*JJT(i).*ones(BEHx,BEHx,1,'uint8'); % same for JJT
    BEH(1:BEHx,2*BEHx+[1:BEHx],[1 3])=zeros(BEHx,BEHx,2,'uint8'); 

    position_sound = 1 + round((i-position(1)).*delta_sound); % for each position (in frame)
    if JJL(i+1)==1, % if there is a lick event
        sound(lick_sound,Fs); % play a lick sound and add it to the sound_band (wav file)
        try, sound_band(position_sound+1:position_sound+length(lick_sound))= ...
            sound_band(position_sound+1:position_sound+length(lick_sound)) + lick_sound; end
    end
    if JJR(i+1)==1, 
        sound(drop_sound,Fs); 
        try, sound_band(position_sound+1:position_sound+length(drop_sound2))= ...
            sound_band(position_sound+1:position_sound+length(drop_sound2)) + drop_sound2; end
    end
    if JJT(i+1)==1, 
        sound(task_sound,Fs); 
        try, sound_band(position_sound+1:position_sound+length(task_sound2))= ...
            sound_band(position_sound+1:position_sound+length(task_sound2)) + task_sound2; end
    end
    
    if varOpen.varOutput.sBe ~= 0 % if the is behavior data, concactenate in x-axis behavior and calcium in III
        III = cat(1,BEH,IMG);
    else % if not, use just calcium
        III = IMG;
    end

    
    
    % --------- 4.2 recording avi file or display -------------------------
    if varOpen.Video.rec ~= 0, % if recording requested  
        if varOpen.Video.rec == 1, % if rec=1 (only image)
            writeVideo(vidObj,III); % store the calcium image in the avi
            disp(['record videoframe ' num2str(i) '/' num2str(position(2))]); % and report where we are in the process 
        elseif varOpen.Video.rec == 2,  % if rec=2 (image+graph) recorded, so snapshot of display
            subplot(121), imshow(III); % in subplot 1: display the image
            title([num2str(myround(i/(varOpen.varOutput.SFr),1)) 's']) % put in title, the time location (in sec)
            subplot(122), % in subplot 2: display the graph
            hold on
            try, delete(h); end % delete the previous time indicator (see next line)
            h = plot([i/(varOpen.varOutput.SFr) i/(varOpen.varOutput.SFr)],[min(seed_info.ca) seed_info.base_beh],'m','LineWidth',2); % display the time indicator
            hold off
            pause(.01)
            F = getframe; % snap show of the figure of the graph         
            RGB = imresize(III,size(F.cdata,1)/size(III,1),'nearest'); % set the y-axis, same size of the behavior image
            RGB = cat(2,RGB,F.cdata); % concatenate in x-axis behavior image and graph
            
            writeVideo(vidObj,RGB); % store the result image in the avi
            disp(['record videoframe (with graph) ' num2str(i)  '/' num2str(position(2))]); % and report where we are in the process 
        end
    else % if no recording, just display
        subplot(121), imshow(III); % in subplot 1: display the image
        title([num2str(myround(i/(varOpen.varOutput.SFr),1)) 's'])  % put in title, the time location (in sec)  
        subplot(122), % in subplot 2: display the graph
        hold on
        try, delete(h); end % see before, time cursor
        h = plot([i/(varOpen.varOutput.SFr) i/(varOpen.varOutput.SFr)],[min(seed_info.ca) seed_info.base_beh],'m','LineWidth',2);  
        hold off
        pause(1/varOpen.Video.fps) % pause at the framerate period
    end
    
end; % end of the time window range
end; % end of the repetition loop


if varOpen.Video.rec ~= 0, 
    close(vidObj); % close the avi file
    disp([varOpen.Video.videofile ' recorded']); 
    filewav = [varOpen.Video.videofile(1:length(varOpen.Video.videofile)-3) 'wav']; % establish wav file name
    audiowrite(filewav,sound_band,Fs); % record the audio file associated
    disp([filewav ' recorded']);

end

