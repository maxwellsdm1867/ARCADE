%A.R.C.A.D.E. stands for ARthur and Chris Analyzing Drosophila Experiments
%#LifeIsLikeAGame


%this version is a proof of concept scirpt, should be modulize later.
close all
clear all
clc


%Mian steps are as follow:
% 0. Specify working directory and load the raw and metadata
% 1. Data parcing and extracting for the csv data file
% 2. Video importing and processing, masking two balls, cropping 
% 3. Dimension reduction using SVD

%% Specify working directory and load the raw and metadata
cd('D:\TuthillLab\data')%go to the data folder, it should contain chris_data.mat and the video you wish to analyze
load('D:\TuthillLab\data\chris_data.mat')%loading the convterd csv table, this contains the metadata and calcuim trace

%% 1. Data parcing and extracting for the csv data file
nov = '20210713_A01_00001_legCam_1.avi'; %name of the video wish to analyze

data = data(2:end,:);%take out the first row, since origianlly it's just text
trail_id =nov(1:18); % trail ID in the data table

%extract the meta data for that video
trail_list = string(data{:,5});% list of trails extracted from the data
trail_target_index = find(trail_list == trail_id);% index of the video we want to analyze
sub_data = data(trail_target_index,:); %extract the subset of the data that is realted to the video trail
calcuim_trace = double(string(sub_data{:,7}));%calcuim traces
frame_number = double(string(sub_data{:,6}));%frame number for alignment
analyze_chek = double(string(sub_data{:,11}));%if zero then don't analyse 
%meta data 
metaData.roi = string(sub_data{1,1});%reigon of interest 
meta.Driver = string(sub_data{1,2});%driver line of the fly 
meta.Ball = string(sub_data{1,3});%on ball or off ball
meta.Animal = string(sub_data{1,4});%animal id 


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
%xbound = [200,570];
xbound = [200,450];
axis image
xlim(ybound)%boundary of x
ylim(xbound)%boundary of y
title('cleaned up image')
colormap('gray')
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
[B,FitInfo] = lasso(A',calcuim_trace,'Alpha',0.5);
Bb= B(:,1)';
U_r = U(:,1:r);
imr = U_r*Bb';
y_hat = (Bb*A)+FitInfo.Intercept(1)*ones(1,length((Bb*A)));
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
cd('D:\TuthillLab\figure\220202')
mkdir(trail_id)
cd(trail_id)
close all
r =10;%take first r pc 
V_r = V_conj(1:r,:);%reduce to r-dimensions
S_r = S(1:r,1:r);%reduce to r-dimensions
A = S_r*V_r;%project to r-dimensions
% lasso

[B,FitInfo] = lasso(A',calcuim_trace,'Alpha',0.5);%elsAtic net
[~,min_idx]= min(FitInfo.MSE);%find the minium lambda values
singular_value_rank = min_idx; 
Bb= B(:,singular_value_rank)';
U_r = U(:,1:r);
imr = U_r*Bb';
y_hat = (Bb*A)+FitInfo.Intercept(singular_value_rank )*ones(1,length((Bb*A)));
%threshold Bb(beta values here)
beta_thrs =0;% 0.003;%thrshold of beta values
Bb_threshold = Bb;
Bb_threshold(abs(Bb)<beta_thrs) = 0;
imr_threshold = U_r*Bb_threshold';
imt_v = 0;
imr_threshold(abs(imr_threshold)<imt_v) = 0;

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
xlabel('difference values')
ylabel('count')
title(['mean of Ca2+ trace' num2str(mean(calcuim_trace))])
saveas(gca,['eph_',num2str(r) '.jpg'])

figure
imagesc(abs(reshape(imr_threshold,size(down_test_frame,1),size(down_test_frame,2))))
axis square
title(['number of the PCs used: ' num2str(r) ' abs(heatmap)'])
colorbar
saveas(gca,['abshm_',num2str(r) '.jpg'])

%model sceletion 
figure
plot(FitInfo.Lambda, FitInfo.MSE)
xlabel('Lambda')
ylabel('MSE')
saveas(gca,['models_',num2str(r) '.jpg'])



figure
sk = histogram(Bb);
title('Beta values distrubtion')
saveas(gca,['bvd_' num2str(r) '.jpg'])


figure
imagesc(douball_mask.*double(ref_frame_m));
ybound = [120, 746];
%xbound = [200,570];
xbound = [200,450];
axis square
xlim(ybound)%boundary of x
ylim(xbound)%boundary of y
title('cleaned up image')
colormap('gray')
saveas(gca,['og_' num2str(r) '.jpg'])


% figure
% imagesc(*abs(reshape(imr_threshold,size(down_test_frame,1),size(down_test_frame,2))))
% axis square
% title(['number of the PCs used: ' num2str(r) ' abs(heatmap)'])
% colorbar
% colormap hot
% saveas(gca,['abshm_',num2str(r) '.jpg'])

% %% error scaling
% in_range = 1:5:100;
% error_vec = error_scaling(in_range,calcuim_trace,A,V_conj,S_r);
% 
% 
% 
% V_r = V_conj(1:r,:);%reduce to r-dimensions
% S_r = S(1:r,1:r);%reduce to r-dimensions
% A = S_r*V_r;%project to r-dimensions
% % lasso
% singular_value_rank = 1; 
% [B,FitInfo] = lasso(A',calcuim_trace);%lasso
% Bb= B(:,singular_value_rank)';
% 
