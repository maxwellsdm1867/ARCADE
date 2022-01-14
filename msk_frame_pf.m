ref_frame = test_frame(:,:,1);
figure
imagesc(ref_frame)
BW1 = edge(ref_frame,'sobel');
imshow(BW1)

% Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = 632;
imageSizeY = 832;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX1 =316;
centerY1 = -500;
radius1 = 650;
circlePixels = (rowsInImage - centerY1).^2 ...
    + (columnsInImage - centerX1).^2 <= radius1.^2;
% circlePixels is a 2D "logical" array.
% Now, display it.
figure
image(circlePixels2) ;
colormap([0 0 0; 1 1 1]);
title('Binary image of a circle');


%image making
ref_frame = test_frame(:,:,1);%reference frame
ref_sizeX = size(ref_frame,1);%length 
ref_sizeY = size(ref_frame,2);%width
imageSizeX = ref_sizeY;
imageSizeY = ref_sizeX;
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX1 =416;
centerY1 = -660;
radius1 = 851;
circlePixels1 = (rowsInImage - centerY1).^2 ...
    + (columnsInImage - centerX1).^2 <= radius1.^2;

centerX2 =280;
centerY2 = 1170;
radius2 = 740;
circlePixels2 = (rowsInImage - centerY2).^2 ...
    + (columnsInImage - centerX2).^2 <= radius2.^2;
mask1 = double(~(circlePixels1|circlePixels2));
figure
imagesc(mask1.*double(ref_frame))
colormap(gray)

ball1= [416,-660,857];
ball2= [280, 1170,740];
rslt =fly_O_msk(ball1,ball2,ref_frame);
figure
imagesc(rslt.*double(ref_frame))