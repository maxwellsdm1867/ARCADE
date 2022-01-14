function rslt =fly_O_msk(ball1,ball2,ref_frame)
%ball1 = [centerX1,centerY1,radius1]
%ball2 = [centerX2,centerY2,radius2]
%ref_frame is for determine size

ref_sizeX = size(ref_frame,1);%length 
ref_sizeY = size(ref_frame,2);%width
imageSizeX = ref_sizeY;
imageSizeY = ref_sizeX;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX1 =ball1(1);
centerY1 = ball1(2);
radius1 = ball1(3);
circlePixels1 = (rowsInImage - centerY1).^2 ...
    + (columnsInImage - centerX1).^2 <= radius1.^2;

centerX2 =ball2(1);
centerY2 = ball2(2);
radius2 = ball2(3);
circlePixels2 = (rowsInImage - centerY2).^2 ...
    + (columnsInImage - centerX2).^2 <= radius2.^2;
mask1 = double(~(circlePixels1|circlePixels2));
rslt = mask1;

end

