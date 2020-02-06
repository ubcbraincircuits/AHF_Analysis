function I2 = OIA_tempobin(I,r)
% function I = OIA_tempobin(I,r)
% resize [X,Y,Z,Z1,Z2] in [X,Y,Z*r,Z1,Z2]
if r ~= 1
    for i = 1:size(I,4)
        for j = 1:size(I,5)
            I2_temp = reshape(single(I(:,:,:,i,j)),size(I,1)*size(I,2),size(I,3));
            I2_temp = imresize(I2_temp,[size(I2_temp,1) size(I2_temp,2)*r],'nearest');
            I2(:,:,:,i,j) = reshape(I2_temp,size(I,1),size(I,2),size(I2_temp,2));
        end
    end
else I2 = I; end;