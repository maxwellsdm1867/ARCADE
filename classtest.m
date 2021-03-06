close all
clear all
clc

figurePath = 'D:\TuthillLab\figure\220302';

pass_condition_name= {'roi' , 'ball','driver','animal'};
pass_condition_value = {'medial','1','JR252','11'};
passCondition(1,1:size(pass_condition_name,2)) = pass_condition_name;
passCondition(2,1:size(pass_condition_value ,2)) = pass_condition_value;

analyse = arcade(figurePath);
[data, metaVariableNames,epochList ]= analyse.eopchListGen(false);
targetID = analyse.dataSelector(epochList,passCondition);
[missFileName,checkedID,trailList] = analyse.epochCheck(epochList,targetID);
reduce_list = trailList(1:5,:);
[gluedDataMatrix,gluedMetaData,douball_mask,ref_frame_m,down_test_frame]= videoReaderGlue(reduce_list,data);
[V_r,S_r,A,U] =  analyse.videoPC(25,gluedDataMatrix);
targertTrace = gluedMetaData.calcuim_trace;%-gluedMetaData.prediction;
[B,FitInfo,imr,y_hat]= analyse.lassoReconstruct(U,A,targertTrace,25);

figure
subplot(2,1,1)
imagesc(douball_mask.*double(ref_frame_m));
ybound = [120, 746];
%xbound = [200,570];
xbound = [200,450];
colorbar
axis image
xlim(ybound)%boundary of x
ylim(xbound)%boundary of y
title('reference image')
subplot(2,1,2)
imagesc(abs(reshape(imr,size(down_test_frame,1),size(down_test_frame,2))))
meta_text =string(strjoin(['roi=' metaData.roi 'ball=' metaData. ball 'animal=' metaData.animal]));
text_dim = [0.3 0.25 0.5 .3 ];
title(['number of the PCs used: ' num2str(r) ' abs(heatmap)'])
axis image
colorbar
colormap('default')
annotation('textbox',text_dim,'String',meta_text,'FitBoxToText','on');%text box


