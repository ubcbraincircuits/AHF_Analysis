function roi = OIA_roidraw(R) 

imagesc(R,[0 max(R(:))]); colormap gray
title('DRAW THE ROI, CLICK OUTSIDE TO FINISH')
i = 1; xi = 2; yi = 2; x = 0;
while ((xi>1) & (yi>1) & (xi<size(R,2)) & (yi<size(R,1)))  
    [xi,yi] = ginput(1);

    if ((xi>1) & (yi>1) & (xi<size(R,2)) & (yi<size(R,1))) 
        x(i) = round(xi);
        y(i) = round(yi);
        if length(x) >= 2
            line([x(i-1) x(i)],[y(i-1) y(i)],'Color',[1 0 0],'LineWidth',1);
        end;
        i = i + 1;
    end;
end;


if length(x) <3
    roi = zeros(size(R));
else
    x(i) = x(1);
    y(i) = y(1);
    roi = roipoly(R,x,y);
end;
