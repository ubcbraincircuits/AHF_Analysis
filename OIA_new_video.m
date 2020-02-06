function [Iout,CC] = OIA_new_video(I,v)
% -------------------------------------------------------------------------
% Matthieu Vanni - Sept 30, 2013
% function [Iout,CC] = OIA_new_video(I,v)
% Toolbox to display and save video
% I = matrix input (X,Y,F,C), up to 3 conditions "C"
% Iout, CC = output matrix
%           Iout: values normalized
%           CC: colorframes (X,Y,F,3)
% v = structure of variable (CAN BE EMPTY):
%     v.videofile: filenameof the video file
%                   default: 'E:\test[TIME].avi'; 
%     v.delta: number of frame for each bloc (temporal binning)
%                   default: 1
%     v.f: gaussian filter (sigma)
%                   default: 0
%     v.level: % of anatomic image when background
%                   default: .5
%     v.param: 
%           when p=1 or 2: prestim period (in frame)
%           when p=3: threshold (%) to show (p =3) 
%                   default: size(I,3)
%     v.mv: mode of display/saving
%           =1 if stored in a avi
%           =2 if stored in a mat 
%           =3 if just play the video in loop
%           =4 if just the output
%                   default: 3
%     v.Iref: background
%           =0 if no surimposition
%           =1 if surimposition on the mean
%           =matrix[X Y] : surimposition on a matrix (can be 3D!)
%                   default: 0
%     v.clipping: 
%           =0 : scale min to max
%           >0 : scale mean +/- clipping*sd
%           [a b] : 
%               a [0 to 1] : cutdown = dont presente a(%) of the bottom
%               b [>0]: clipping
%           [a b c d] or [a b c d e f] : idem but for 2 or 3 colors 
%                   default: 0
%     v.colormap: lookup table
%                   default: gray 
%     v.p: processing
%           =0 if nothing
%           =1 if DF/F
%           =2 if DF alone
%           =3 if gradient
%                   default: 0
%     v.fps: frame/s
%                   default: 30 
%     v.size: spatial binning
%                   default: 1
%     v.masque: black mask
%                   default: 0
%     v.SF: sampling frequency (to present the time)
%                   default: 1
%     v.time: position of the time displayed
%                   default: 0 (not displayed)
%     v.L: witdh of the window (to present the scale bar)
%                   default: 1
%     v.bar: position of the bar displayed
%                   default: 0 (not displayed)
%     v.lut: position of the LUT table
%                   default: 0 (not displayed)
%     v.point: vector [frame 2] or [X Y] of position of a point
%                   default: 0 (not displayed)
%     v.colpoint: color of the point [R G B]
%                   default: black [0 0 0]

% ------------------------------- SET DEFAULT VALUES ----------------------

I = single(I);

try, length(v); catch, v.empty = 1; end;

if isfield(v,'clipping') == 1, clipping = v.clipping;
else, clipping = 0; end;
if isfield(v,'param') == 1, param = v.param;
else, param = size(I,3); end;
if isfield(v,'delta') == 1, delta = v.delta;
else, delta = 1; end;
if isfield(v,'f') == 1, f = v.f;
else, f = 0; end;
if isfield(v,'p') == 1, p = v.p;
else, p = 0; end;
if isfield(v,'Iref') == 1, Iref = v.Iref;
else, Iref = 0; end;
if isfield(v,'level') == 1, level = v.level;
else, level = .5; end;
if isfield(v,'fps') == 1, fps = v.fps;
else, fps = 30; end;
if isfield(v,'videofile') == 1, videofile = v.videofile;
else, videofile = ['E:\test' num2str(round(100000000000*now)) '.avi']; end;
if isfield(v,'colormap') == 1, colormap_lut = v.colormap;
else, colormap_lut = gray; end;
if isfield(v,'size') == 1, size_video = v.size;
else, size_video = 1; end;
if isfield(v,'mv') == 1, mv = v.mv;
else, mv = 3; end;
if isfield(v,'masque') == 1, masque = v.masque;
else, masque = 0; end;
if isfield(v,'SF') == 1, SF = v.SF;
else, SF = 1; end;
if isfield(v,'time') == 1, timepos = v.time;
else, timepos = 0; end;
if isfield(v,'L') == 1, L = v.L;
else, L = 1; end;
if isfield(v,'bar') == 1, barpos = v.bar;
else, barpos = 0; end;
if isfield(v,'lut') == 1, lutpos = v.lut;
else, lutpos = 0; end;
if isfield(v,'point') == 1, point = round(v.point); ispoint = 1; 
else, ispoint = 0; end;
if isfield(v,'colpoint') == 1, colpoint = v.colpoint;  
else, colpoint = [0 0 0]; end;



% ------------------------------- SPATIAL RESIZING/BINNING ----------------
size_video_on = 0;
if length(size_video) == 1
    if size_video ~= 1
        size_video_on = 1;
    end;
end
if length(size_video) > 1
    size_video_on = 1;
end
if size_video_on == 1   
    disp('  > resizing...')
    for i = 1:size(I,3)
        for j = 1:size(I,4)
            I2(:,:,i,j) = imresize(I(:,:,i,j),size_video);
        end;
    end
    I = I2; clear I2;
    if size(masque,1) ~= size(I,1)
        masque = imresize(masque,size_video);
    end;
end;
if length(masque) == 1
    masque = ones(size(I,1),size(I,2));
end

% ------------------------------- SET THE BACKGROUND FIELD ----------------
disp('  > set backgroud field...')
incru = 0;
if length(Iref)== 1
    if Iref == 1
        Iref = mean(I(:,:,:,1),3);
        Iref = Iref - min(Iref(:));
        Iref = level*Iref / max(Iref(:));
        incru = 1;
    end;   
else
    for i = 1:size(Iref,3)
        tampon = Iref(:,:,i);
        tampon = tampon - min(tampon(:));
        Iref(:,:,i) = level*tampon / max(tampon(:));
    end;
    incru = 1;
end;

% ------------------------------- TEMPORAL BINNING ------------------------
if delta ~= 1
    disp('  > temporal binning...')
    param = param / delta; % change the "prestim" value
    for j = 1:size(I,4)
        III = I(:,:,:,j);
        J = reshape(III,size(I,1)*size(I,2),size(I,3),1); 
        J = imresize(J,[size(J,1),round(size(J,2)/delta)]);
        JJ(:,:,:,j) = reshape(J,size(I,1),size(I,2),size(J,2));
    end
    I = JJ; clear III JJ
end
tempotag = zeros(size(I,3),1);
tempotag(param:param+3)=ones(4,1);
% ------------------------------- SPATIAL FILTERING -----------------------
if f ~= 0
    disp('  > spatial filtering...')
    f = fspecial('gaussian',5*f,f);
    for j = 1:size(I,4)
        for i = 1:size(I,3)
            I(:,:,i,j) = conv2(I(:,:,i,j),f,'same')./conv2(ones(size(I(:,:,i,j))),f,'same');
        end;
    end
end;

% ------------------------------- PROCESS THE RESPONSES -------------------
disp('  > processing the response...')
clear D
if p == 1 % dans le cas de p == 1, I reste le meme (In the case that p==1, I remains the same)
    for j = 1:size(I,4) 
        Dm = mean(I(:,:,1:param,j),3);
        for i = 1:size(I,3)
            D(:,:,i,j) = 100*(I(:,:,i,j)-Dm)./Dm;  
        end; 
    end
    I = D;
elseif p == 2
    for j = 1:size(I,4) 
        Dm = mean(I(:,:,1:param,j),3);
        for i = 1:size(I,3)
            D(:,:,i,j) = (I(:,:,i,j)-Dm);  
        end; 
    end
    I = D;
elseif p == 3
    for j = 1:size(I,4)  
        for i = 1:size(I,3)-1
            D(:,:,i,j) = I(:,:,i,j)-I(:,:,i+1,j);  
        end;
    end;
    I = abs(D);
end;

% ------------------------------- CLIPPING --------------------------------
disp('  > clipping...')
if length(clipping) == 1
    cutdown = 0;
elseif length(clipping) == 2
    cutdown = clipping(1);
    clipping = clipping(2);
elseif length(clipping) == 4
    cutdown(1) = clipping(1);
    clipping_tampon(1) = clipping(2);
    cutdown(2) = clipping(3);
    clipping_tampon(2) = clipping(4);
    clipping = clipping_tampon;
elseif length(clipping) == 6
    cutdown(1) = clipping(1);
    clipping_tampon(1) = clipping(2);
    cutdown(2) = clipping(3);
    clipping_tampon(2) = clipping(4);
    cutdown(3) = clipping(5);
    clipping_tampon(3) = clipping(6);
    clipping = clipping_tampon;
end;
if ((size(I,4)==2) & (length(clipping)==1))
    clipping(2) = clipping(1);
    cutdown(2) = cutdown(1);
    disp('corr')
end;
if ((size(I,4)==3) & (length(clipping)==1))
    clipping(2) = clipping(1);
    cutdown(2) = cutdown(1);
    clipping(3) = clipping(1);
    cutdown(3) = cutdown(1);
end;
for i = 1:size(I,4)
    D = I(:,:,:,i);
    if clipping(i) > 0
        mD = mean(D(:));
        sD = std(D(:));
        mini = mD-clipping(i)*sD;
        maxi = mD+clipping(i)*sD;
        D = D.*((D>=mini)&(D<=maxi)) + mini.*(D<mini) + maxi.*(D>maxi); 
        D(1,1,:)=mini*ones(1,1,size(D,3));
        D(1,2,:)=maxi*ones(1,1,size(D,3));
    end;
       
    minval = min(D(:));
    maxval = max(D(:));
    diff = (maxval - minval).*(1-cutdown); % used for the display
    
    D = D - min(D(:));
    D = D / max(D(:));
    if p == 3, D = D.*(D>param); end;
    D = D.*(D>cutdown(i))-cutdown(i); 
    D = D.*(D>0);   
    I(:,:,:,i) = D / max(D(:)); 
end
Iout = I;
disp(['      (range value displayed: ' num2str(maxval - diff) '-' num2str(maxval) ')'])



% ------------------------------- LUT RANGE DISPLAY -----------------------
if length(lutpos)>1 
    hLUT = round(size(I,2)/2);
    LUT = OIA_color_intern(imresize([64:-1:1]',[hLUT 3]),jet);

    Text1 = sprintf([num2str(maxval)]);
    Text1 = Text1(1:4);
    H1 = vision.TextInserter(Text1);
    H1.Color = 255*[1 1 1];
    H1.FontSize = 10;
    H1.Location = [lutpos(1) lutpos(2)+3];
    H1.Font = 'Arial';

    Text2 = sprintf([num2str(maxval - diff)]);
    Text2 = Text2(1:4);
    H2 = vision.TextInserter(Text2);
    H2.Color = 255*[1 1 1];
    H2.FontSize = 10;
    H2.Location = [lutpos(1)+hLUT-H2.FontSize lutpos(2)+3];
    H2.Font = 'Arial';
 
end;

% -------------------------- CREATE IMAGE AND SAVE THEM IN AVI FILE -------
disp(['  > create video... check location of varOpen.working_folder for the AVI file'])
if mv == 1, mov = VideoWriter(videofile); mov.FrameRate = fps; open(mov); end;
%if mv == 1, mov = avifile(videofile,'Compression','None','fps',fps); end;
if ispoint == 1
    if size(point,1) == 1
        point = imresize(point,[size(I,3) 2]);
    end;
end;

wait_bar = waitbar(0,'build the video...');
for i = 1:size(I,3) % for each frame (to add in the avi file)
    waitbar(i/size(I,3), wait_bar);
    
    %disp(['frame: ' num2str(i) '/' num2str(size(I,3))])
    clear C
    if incru == 0
        if size(I,4) == 1    
            C = OIA_color_intern(I(:,:,i),colormap_lut); 
        elseif size(I,4) == 2
            C(:,:,1) = I(:,:,i,1);
            C(:,:,2) = I(:,:,i,2);
            C(:,:,3) = zeros(size(I(:,:,1,1)));
            C = uint8(255*C);
        elseif size(I,4) == 3
            for k = 1:3, 
                C(:,:,k) = I(:,:,i,k); 
            end;
            C = uint8(255*C);
        end;    
    else
        if size(Iref,3) == 1
            for k = 1:3, C(:,:,k) = Iref; end;
        else
            for k = 1:3, C(:,:,k) = Iref(:,:,i); end;
        end;
        
        q = sum(colormap_lut);
        if q(1)==q(2) % strange way to identify if this is a gray colorcode (red=green=blue)
            for j = 1:size(I,4) 
                C(:,:,j) = C(:,:,j) + (1-level)*I(:,:,i,j);
            end   
        else
            addC = single(OIA_color_intern(I(:,:,i),colormap_lut))/255; 
            for j = 1:3 
                C(:,:,j) = (((1/level)*C(:,:,j)).*(I(:,:,i)==0)) + (addC(:,:,j).*(I(:,:,i)>0));
                
                
                
            end  
        end;
        C = uint8(255*C);
    end;
    
    for j = 1:3
        C(:,:,j) = C(:,:,j).*uint8(masque);
    end;
    
    
    if length(timepos)>1 % tag of time
        temps = round(10*delta*i/SF)/10;
        Text = sprintf([num2str(temps) ' s']);
        H = vision.TextInserter(Text);
        if mean(masque) == 1, H.Color = 255*[0 0 0];
        else, H.Color = 255*[1 1 1]; end
        H.FontSize = 10;
        H.Location = timepos;
        H.Font = 'Arial';
        C = step(H, C);
    end;
    if length(lutpos)>1
        C(lutpos(1)+1:lutpos(1)+hLUT,lutpos(2):lutpos(2)+2,:) = LUT;
        C = step(H1, C);
        C = step(H2, C);
    end;
    
    if length(barpos)>1 % scale bar 
        SB =round(size(C,1)/L);
        C(barpos(1),barpos(2)+1:barpos(2)+SB,:) = 255*ones(1,SB,3,'uint8');
    end
    
    
    % ------------------------------- TAG THE STIMULUS ------------------------
    if ((p == 0) | (p == 1) | (p == 2))
        if tempotag(i)==1    
            C(1:round(size(C,1)/10),1:round(size(C,2)/10),:) = 255*ones(round(size(C,1)/10),round(size(C,2)/10),3,'uint8');
        end
        
    end;

    if ispoint == 1
        C(round(point(i,1)),round(point(i,2)),:) = 255*reshape(colpoint,1,1,3);
    end;
    


    if mv==1
        writeVideo(mov,C);
        %mov = addframe(mov,C); 
    end
    CC(:,:,:,i) = C; 
end;              
close(wait_bar);

if mv == 1 
    close(mov); 
elseif mv == 2
    save(videofile,'Iout','CC')
elseif mv == 3
    for r = 1:1000
        for i = 1:size(CC,4)
            imshow(CC(:,:,:,i))
            title(num2str(i))
            pause(1/fps)
        end;
        pause(1/(fps/5))
    end;
end;

function J = OIA_color_intern(T,LUT_intern)
T = T - min(T(:));
T = round(255 * T / max(T(:)));
LUT_intern = imresize(LUT_intern,[256 3]);
J = 255*ind2rgb(T,LUT_intern);
J = uint8(J);


