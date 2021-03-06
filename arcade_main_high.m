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


%modifications
%%add in new tweaks 20220309
%1. use SVDs
%2. low pass and high pass the calcuim data
%3. walkig/non-walking extraction
%4. piece together for one amimal
%% Specify working directory and load the raw and metadata
cd('D:\TuthillLab\video')%go to the data folder, it should contain chris_data.mat and the video you wish to analyze
load('D:\TuthillLab\video\data')%loading the convterd csv table, this contains the metadata and calcuim trace

%% 1. Data parcing and extracting for the csv data file
nov = '20211103_A01_00020_legCam_1.avi'; %name of the video wish to analyze

data = data(2:end,:);%take out the first row, since origianlly it's just text
trail_id =nov(1:18); % trail ID in the data table

%extract the meta data for that video
trail_list = string(data{:,5});% list of trails extracted from the data
trail_target_index = find(trail_list == trail_id);% index of the video we want to analyze
sub_data = data(trail_target_index,:); %extract the subset of the data that is realted to the video trail
calcuim_trace_pre = double(string(sub_data{:,7}));%calcuim traces
frame_number_pre = double(string(sub_data{:,6}));%frame number for alignment
analyze_chek = double(string(sub_data{:,12}));%if zero then don't analyse
analyze_flag = sum(analyze_chek)/length(analyze_chek);%right percentage
walk_flag = double(string(sub_data{:,14}));
prediction_pre = double(string(sub_data{:,8}));
rest_flag = double(string(sub_data{:,13}));
flag_pole = analyze_chek;%.*(walk_flag==0);%this is and opertaion 
calcuim_trace = calcuim_trace_pre(flag_pole==1);%clean the calcuim trace
prediction_trace  = prediction_pre(flag_pole==1);%clean up the prediction
frame_number = frame_number_pre(flag_pole==1);%clean the frame number

%meta data
metaData.roi = string(sub_data{1,1});%reigon of interest
metaData.driver = string(sub_data{1,2});%driver line of the fly
metaData.ball = string(sub_data{1,3});%on ball or off ball
metaData.animal = string(sub_data{1,4});%animal id
metaData.id = trail_id; % trail id number
metaData%display meta data

%% calcuim trace processing
trace_difference = calcuim_trace-prediction_trace;

%% Video importing and processing
cd('C:\Users\USER\Downloads\Data')
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
%sanity check figure
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

%read the video frame by frame
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

%take 25pcs
k_components = 25;
tic
[U,S,V] = svds(video_matrix,k_components);%SVD of data matrix
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

%% 
rect_difference = trace_difference;
for i=1:length(rect_difference)
if rect_difference(i)<0
    rect_difference(i) =0;
end
end
target_calcuim = calcuim_trace;


%debugging
cd('D:\TuthillLab\figure\220309')
mkdir(trail_id)
cd(trail_id)
close all
r =k_components;%take first r pc
V_r = V_conj(1:r,:);%reduce to r-dimensions
S_r = S(1:r,1:r);%reduce to r-dimensions
A = S_r*V_r;%project to r-dimensions
% lasso
[B,FitInfo] = lasso(A',target_calcuim,'Alpha',0.5);%elsAtic net
[~,min_idx]= min(FitInfo.MSE);%find the minium lambd
singular_value_rank = min_idx; %use the correspond singular value
Bb= B(:,singular_value_rank)';%get n-th row out


U_r = U(:,1:r);%reduce U to r by b
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

%add in textbox to display the meta data%TO DO

figure
hold on
plot(target_calcuim)
plot(y_hat')
hold off
legend('Ca2+ trace','model prediction')
title(['number of the PCs used: ',num2str(r)])
saveas(gca,['MP_' num2str(r) '_h.jpg'])

figure
err_pre = y_hat-calcuim_trace';
kk = histogram(err_pre);
xlabel('difference values')
ylabel('count')
title(['mean of Ca2+ trace' num2str(mean(calcuim_trace))])
saveas(gca,['eph_',num2str(r) '.jpg'])

figure
hold on 
plot(calcuim_trace)
plot(prediction_trace)
plot(trace_difference)
hold off
title('calcuim check')
legend('calcuim trace','prediction trace','trace difference')
saveas(gca,['cc_',num2str(r) '.jpg'])
%model sceletion
% figure
% hold on
% plot(FitInfo.MSE,sparsity)
% scatter(FitInfo.MSE(singular_value_rank ),sparsity(singular_value_rank ))
% ylabel('%of zeros in beta')
% xlabel('MSE')
% hold off
% saveas(gca,['models_',num2str(r) '.jpg'])



figure
sk = histogram(Bb);
title('Beta values distrubtion')
saveas(gca,['bvd_' num2str(r) '.jpg'])


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
imagesc(abs(reshape(imr_threshold,size(down_test_frame,1),size(down_test_frame,2))))
meta_text =string(strjoin(['roi=' metaData.roi 'ball=' metaData. ball 'animal=' metaData.animal]));
text_dim = [0.3 0.25 0.5 .3 ];
title(['number of the PCs used: ' num2str(r) ' abs(heatmap)'])
axis image
colorbar
colormap('default')
annotation('textbox',text_dim,'String',meta_text,'FitBoxToText','on');%text box
saveas(gca,['og_' num2str(r) '_h.jpg'])

