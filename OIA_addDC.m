function If=OIA_addDC(If,DC)
% function If=OIA_addDC(If,DC)
% add DC matric to each frame
%disp('  > re-add DC component')
for i = 1:size(If,3)
    If(:,:,i) = single(If(:,:,i)) + single(DC);
end;
