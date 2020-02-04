function FARM_singletrial(varOpen) 

qtype_list = varOpen.singletrial.qtype_list; % determine qtype 

qtype_files= {};
for i= 1:length(qtype_list) 
    qtype_files{i} = ['qtype' num2str(qtype_list(i)) '.mat'], % extract the averaging data from a specific qtype 
end 

qtype_pic_files= {} ;
for i = 1:length(qtype_list) 
    qtype_pic_files{i} = ['mouse' varOpen.mouse '_qtype' num2str(qtype_list(i)) '_trials' num2str(min(varOpen.singletrial.trial_interval)) '-' num2str(max(varOpen.singletrial.trial_interval)) '.pdf']; % save the data as pdf 
end

for m = 1:length(qtype_list) 
    
load(qtype_files{m},'I2','I2_be') % first input is the name of the qtype file 

    %keep only the GCaMP channel
I=I2(:,:,:,1,:); % viewing the fluorescence data(channel = 1); can change it to 2nd channel for reflectance
clear I2
I=squeeze(I);

%        I2
% keep only the needed events (qtype1)
% 1st input indicate x dimension of frame 
% 2nd input indicates y dimension
% 3rd input indicates the time 
% 4th represents fluorescence(1) and reflectance(2) channels 
% 5th represents the number of trials (probably will not put all trials on there because it's hard to see )

% I=I(:,:,:,1:58);
if length(varOpen.singletrial.trial_interval) == 0 
    varOpen.singletrial.trial_interval = [1:size(I2_be,3)] ; 
else 
end 

I=I(:,:,:,varOpen.singletrial.trial_interval); 

%     I2_be 
% 1st input shows the number of ROIs selected 
% 2nd shows the time 
% 3rd shows the number of trials 
% same for behaviour cam

%      I4 
% 1st input shows x dimension of frame 
% 2nd shows the y dimension of frame 
% 3rd shows the time 

be=I2_be(:,:,varOpen.singletrial.trial_interval);
clear I2_be
%keep only the roi of interest
be_all=be;
be=squeeze(be(varOpen.singletrial.roi,:,:)); 

%added by Jeff to allow user to use [] input to signify the whole time rage
if isempty(varOpen.singletrial.time)    
varOpen.singletrial.time=[1 size(I,3)];
end

%more temporal binning
I=I(:,:,varOpen.singletrial.time(1):varOpen.singletrial.time(2),:); % taking into account frames user wants to view 
% determining size(I,4) gives you the number of trials you are evaluating 
% I=reshape(I,64,64,varOpen.singletrial.temp_bin,(size(I,3)/varOpen.singletrial.temp_bin),size(I,4));
I=reshape(I,64,64,((length(varOpen.singletrial.time(1):varOpen.singletrial.time(2)))/varOpen.singletrial.temp_bin),varOpen.singletrial.temp_bin,size(I,4)); % Taking the time and splitting it into two dimensions (1:80 -> 2,40) 
I=mean(I,4); % take the mean of every column to generate the number of images specified in the 4th dimension of I stated above, generating I(:,:,1,n of images, n of trials) 
%teri had a 3 above (ln62) which was causing the binning to be done on the wrong
%index.
I=squeeze(I); % take out the 3rd dimension
I_ave= mean(I,4);
I_ave= squeeze(I_ave);
% size of I = (xdim, ydim, n of images, n of trials) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAKING THE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure;
fig1.Renderer = 'Painter'; 
% make subplots for each trial
ntrials=size(be,2); 
nplots=0; % number of plots for the brain montages on the left hand side 
clims=[0 10]; % colour limit 
for i=1:ntrials
    %first plot the montage
    nplots=nplots+1; % placing frame of nth trial to (n-1)th trial 
    subtightplot(ntrials,2,nplots,[.001,.001]) % giving it little to no space between the frames 
    to_plot=I(:,:,:,i);% taking the image dimensions and time data into account for a specific trial 
    % changing the y dimension every time, xdim stays the same but to put the plots next to each other, combine the previous plot with the
    % current to create one new plot with larger ydim. 
    % final size of montage = xdim by ydim*n images
    to_plot=reshape(to_plot,(size(to_plot,1)),size(to_plot,2)*size(to_plot,3)); 
    imagesc(to_plot,clims);
%     xlim([0 640])
    set(gca,'XTick',[], 'YTick', [])% remove x and y axis 
    axis image;
    colormap jet
    nplots=nplots+1; % this plots the roi curve next to the montage by adding 1 to nplots 
    subtightplot(ntrials,2,nplots,[.0001,.0001])
    plot(be(:,i),'r'); % size of be = time x ntrials
    xlim([0 40]) % (time) range of movement you want to see (window of data you want to view) 
    set(gca,'XTick',[], 'YTick', []) 
    switch qtype_list(m) % if positive, this is one of the text tage
           case 2, qstring = 'GO=2'; % look for 'GO=2' in text
           case 3, qstring = 'GO=-4';
           case 4, qstring = 'GO=1';
           case 5, qstring = 'GO=-1';
           case 6, qstring = 'GO=-3';
           otherwise, qstring = 'noanychancethisstringcanexist';
       end
    
     p= mtit(['Mouse ' varOpen.mouse ' ' qstring ' For Trials ' num2str(min(varOpen.singletrial.trial_interval)) '-' num2str(max(varOpen.singletrial.trial_interval))]) ; % title 
end
qtype_single_figure= [varOpen.singletrial.qtypeoutfile qtype_pic_files{m} ]; disp(qtype_single_figure)
saveas(gcf,qtype_single_figure)

tmpbx=varOpen.varOutput.bx;tmpby=varOpen.varOutput.by;
save([varOpen.working_folder 'single_trial_data.mat'],'I','be','tmpbx','tmpby','be_all');

figure;
nplots = 0;
for i = 1 
    nplots=nplots+1; % placing frame of nth trial to (n-1)th trial 
    subtightplot(ntrials,2,nplots,[.001,.001]) % giving it little to no space between the frames 
    to_plot=I_ave(:,:,:); % taking the image dimensions and time data into account for a specific trial 
    % changing the y dimension every time, xdim stays the same but to put the plots next to each other, combine the previous plot with the
    % current to create one new plot with larger ydim. 
    % final size of montage = xdim by ydim*n images
    to_plot=reshape(to_plot,(size(to_plot,1)),size(to_plot,2)*size(to_plot,3)); 
    imagesc(to_plot,clims);
%     xlim([0 640])
    set(gca,'XTick',[], 'YTick', [])% remove x and y axis 
    axis image;
    colormap jet
    nplots=nplots+1; % this plots the roi curve next to the montage by adding 1 to nplots 
    subtightplot(ntrials,2,nplots,[.001,.001])
    plot(be(:,1),'r'); % size of be = time x ntrials
    xlim([0 40]) % (time) range of movement you want to see (window of data you want to view) 
    set(gca,'XTick',[], 'YTick', []) 
    
        switch qtype_list(m) % if positive, this is one of the text tage
           case 2, qstring = 'GO=2'; % look for 'GO=2' in text
           case 3, qstring = 'GO=-4';
           case 4, qstring = 'GO=1';
           case 5, qstring = 'GO=-1';
           case 6, qstring = 'GO=-3';
           otherwise, qstring = 'noanychancethisstringcanexist';
       end
    
     p= mtit(['Average of Mouse ' varOpen.mouse ' ' qstring ' For Trials ' num2str(min(varOpen.singletrial.trial_interval)) '-' num2str(max(varOpen.singletrial.trial_interval))]) ;
end
qtype_single_figure= [varOpen.singletrial.qtypeoutfile 'average' qtype_pic_files{m} ];
saveas(gcf,qtype_single_figure)
clear nplots ntrials to_plot i 
end

end

