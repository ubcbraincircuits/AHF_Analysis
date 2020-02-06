function I=OIA_multiDC(If,DC)
% function If=OIA_multiDC(If,DC)
% multiply DC matric to each frame
%disp('  > re-add DC component')
dimen = length(size(If));
if dimen==3
    for i = 1:size(If,3)
        I(:,:,i) = single(If(:,:,i)) .* single(DC);
    end
elseif dimen==4
    for i = 1:size(If,3)
        for j = 1:size(If,4)
            I(:,:,i,j) = single(If(:,:,i,j)) .* single(DC);
        end
    end;
elseif dimen==5
    for i = 1:size(If,3)
        for j = 1:size(If,4)
            for k = 1:size(If,5)
                I(:,:,i,j,k) = single(If(:,:,i,j,k)) .* single(DC);
            end
        end
    end;
end
