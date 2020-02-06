function J = OIA_color(T,LUT)
% ------------------------------------
% function J = OIA_color(I,LUT) ou
% function J = OIA_color(I)
% Matthieu Vanni - 28 aout 2004
% Transforme une image [x y 1] en une image couleur
% 24 bits [x y 3] (3x8 bits)
% I = image [x y 1]
% J = image [x y 3]
% facultatif : LUT (par defaut = gray)
% NB : autre colormap =
% RAY, HOT, COOL, BONE, COPPER, 
% PINK, FLAG, PRISM, JET, GRAY, HSV
% ------------------------------------
if nargin == 1
    LUT = gray; 
end;

Nanloc = isnan(T);
Nanmini = nanmin(T(:));
for i = 1:size(Nanloc,1)
    for j = 1:size(Nanloc,2)
        if Nanloc(i,j)==1, T(i,j)=Nanmini; end
    end
end
T = T - min(T(:));
T = round(255 * T / max(T(:)));
LUT = imresize(LUT,[256 3]);
J = 255*ind2rgb(T,LUT);
for i = 1:3,
    J(:,:,i) = J(:,:,i).*~Nanloc;
end
J = uint8(J);
