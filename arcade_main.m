%A.R.C.A.D.E. stands for ARthur and Chris Analyzing Drosophila Experiments
%#LifeIsLikeAGame


%this version is a proof of concept scirpt, should be modulize later.
close all
clear all
clc


%Mian steps are as follow:
% 0. Specify working directory and load the raw and metadata
% 1. Data parcing and extracting for the csv data file
% 2. Video importing and processing, masking two balls
% 3. Dimension reduction using SVD

%% Specify working directory and load the raw and metadata
cd('D:\TuthillLab\data')%go to the data folder, it should contain chris_data.mat and the video you wish to analyze
load('D:\TuthillLab\data\chris_data.mat')%loading the convterd csv table, this contains the metadata and calcuim trace

%% 1. Data parcing and extracting for the csv data file
nov = '20210716_A01_00002_legCam_1.avi'; %name of the video wish to analyze

data = data(2:end,:);%take out the first row, since origianlly it's just text
trail_id =nov(1:18); % trail ID in the data table

%extract the meta data for that video
trail_list = string(data{:,5});% list of trails extracted from the data
trail_target_index = find(trail_list == trail_id);% index of the video we want to analyze
sub_data = data(trail_target_index,:); %extract the subset of the data that is realted to the video trail
calcuim_trace = double(string(sub_data{:,7}));%calcuim traces
frame_number = double(string(sub_data{:,6}));%frame number for alignment

%% Video importing and processing
vod_read = VideoReader(nov); %read the correspondence video
frame_start = frame_number(1);%start frame
frame_end = frame_number(end);%end frame
%compress the video
compress_factor = 0.5;
test_frame= read(vod_read,100);
%masking 
ball1= [416,-660,857];%upper ball [centerX1,centerY1,radius1]
ball2= [280, 1290,860];%lower ball [centerX1,centerY1,radius1]
ref_frame_m = test_frame(:,:,1);
douball_mask =fly_O_msk(ball1,ball2,ref_frame_m);% mask image
figure
imagesc(douball_mask.*double(ref_frame_m));
ybound = [120, 746];
xbound = [200,570];
xlim(ybound)%boundary of x
ylim(xbound)%boundary of y
title('cleaned up image')
im_crop = ref_frame_m(xbound(1):xbound(2),ybound(1):ybound(2));
down_test_frame = imresize(im_crop, compress_factor);
video_matrix = zeros(size(down_test_frame,1)*size(down_test_frame,2),length(frame_number));% data matrix, each column represent one frame
parfor i = 1:length(frame_number) %use parallel worker for faster computing 
    frame_idx = frame_number(i);
    temp_frame = read(vod_read,frame_idx);
    temp_frame1 = double(temp_frame(:,:,1)).*douball_mask;%reduce to 2 dimension and mask;
    %add down size here
    temp_frame12 = temp_frame1(xbound(1):xbound(2),ybound(1):ybound(2));
    temp_frame2 = imresize(temp_frame12, compress_factor);%video compressing
    video_matrix(:,i)= reshape( temp_frame2,[],1);% this shold obtain the frames in their order
end

%sainity check
% figure
% for k = 500:15000
%     pause(0.0005)
%  imagesc(reshape(video_matrix(:,k),size(down_test_frame,1), size(down_test_frame,2)))
% end
% 

%% Dimension reduction using SVD

%for 25% compression it takes around 550 seconds 
tic
[U,S,V] = svd(video_matrix);%SVD of data matrix
toc

%try to accrealte this later 
% tic
% f= parfeval(@svd,3,video_matrix);
% [U1,S1,V1]=fetchOutputs(f);
% toc


sv = diag(S);
figure
stairs((cumsum(sv)/sum(sv)).*(cumsum(sv)/sum(sv)))
title('Explained varience')
V_conj = V';
r = 100;%take first r pcs 
V_r = V_conj(1:r,:);%reduce to r-dimensions
%V_r = V_conj(:,1:r);%reduce to r-dimensions

S_r = S(1:r,1:r);%reduce to r-dimensions
A = S_r*V_r;%project to r-dimensions
%% lasso
[B,FitInfo] = lasso(A',calcuim_trace);
Bb= B(:,1)';
U_r = U(:,1:r);
imr = U_r*Bb';
y_hat = (Bb*A)+FitInfo.Intercept(1)*ones(1,14527);
figure
hold on 
plot(calcuim_trace)
plot(y_hat')
hold off 
legend('Ca2+ trace','model prediction')
title(['number of the pc used ',num2str(r)])
figure
err_pre = y_hat-calcuim_trace;
kk = histogram(err_pre);
title(['mean of Ca2+ trace' num2str(mean(calcuim_trace))])
figure
imagesc(down_test_frame);
colormap(gray)
figure
imagesc(abs(reshape(imr,size(down_test_frame,1),size(down_test_frame,2))))



%%
%debugging
cd('D:\TuthillLab\figure\220119')
close all
r = 100;%take first r pc 
V_r = V_conj(1:r,:);%reduce to r-dimensions
S_r = S(1:r,1:r);%reduce to r-dimensions
A = S_r*V_r;%project to r-dimensions
% lasso
singular_value_rank = 1; 
[B,FitInfo] = lasso(A',calcuim_trace);%lasso
Bb= B(:,singular_value_rank)';
U_r = U(:,1:r);
imr = U_r*Bb';
y_hat = (Bb*A)+FitInfo.Intercept(1)*ones(1,14527);

%plotting
figure
hold on 
plot(calcuim_trace)
plot(y_hat')
hold off 
legend('Ca2+ trace','model prediction')
title(['number of the PCs used: ',num2str(r)])
saveas(gca,['MP_' num2str(r) '.jpg'])

figure
err_pre = y_hat-calcuim_trace';
kk = histogram(err_pre);
title(['mean of Ca2+ trace' num2str(mean(calcuim_trace))])
saveas(gca,['eph_',num2str(r) '.jpg'])

figure
imagesc(abs(reshape(imr,size(down_test_frame,1),size(down_test_frame,2))))
title(['number of the PCs used: ' num2str(r) ' abs(heatmap)'])
colorbar
saveas(gca,['abshm_',num2str(r) '.jpg'])


figure
imagesc((reshape(imr,size(down_test_frame,1),size(down_test_frame,2))))
colorbar
title(['number of the PCs used: ' num2str(r) ' heatmap'])
saveas(gca,['hm_' num2str(r) '.jpg'])
