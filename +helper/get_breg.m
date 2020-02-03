function [bx,by] = get_breg(img)
% Get Bregma Coordinates
%
% Input:  img      (Image data. Use fluorescence)
% Output: [bx,by]  (x and y location of Bregma in pixels)
% 
% Usage:  [bx,by] = get_breg(img);

figure, imagesc(img(:,:,1)), colormap gray 
title('Click on Bregma')

[bx,by] = ginput(1);

%close gcf
