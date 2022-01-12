%A.R.C.A.D.E. stands for ARthur and Chris Analyzing Drosophila Experiments
%#LifeIsLikeAGame


%this version is a proof of concept scirpt, should be modulize later.
close all
clear all
clc


%Mian steps are as follow:
% 0. Specify working directory and load the raw and metadata
% 1. Data parcing and extracting for the csv data file
% 2. Video importing and processing
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

compress_factor = 0.25;
test_frame= read(vod_read,1);
down_test_frame = imresize(test_frame(:,:,1), compress_factor);


video_matrix = zeros(size(down_test_frame,1)*size(down_test_frame,2),length(frame_number));% data matrix, each column represent one frame
parfor i = 1:length(frame_number) %use parallel worker for faster computing 
    frame_idx = frame_number(i);
    temp_frame = read(vod_read,frame_idx);
    temp_frame1 = temp_frame(:,:,1);%reduce to 2 dimension;
    temp_frame2 = imresize(temp_frame1, compress_factor);%video compressing
    video_matrix(:,i)= reshape( temp_frame2,[],1);% this shold obtain the frames in their order
end

%sainity check
% figure
% for k = 500:15000
%     pause(0.0005)
%  imagesc(reshape(video_matrix(:,k),frame_size(1), frame_size(2)))
% end


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
stairs(cumsum(sv)/sum(sv))
V_conj = V';
r = 5000;%take first 5000 pc 
V_r = V(1:r,:);%reduce to r-dimensions
S_r = S(1:r,1:r);%reduce to r-dimensions
A = S_r*V_r;%project to r-dimensions

%% lasso
B = lasso(A',calcuim_trace');
Bb= B(:,1)';
k = U(:,1:r);
imr = k*Bb';
figure
imagesc(down_test_frame)
colormap(gray)
hold on
imagesc(reshape(imr,size(down_test_frame,1),size(down_test_frame,2)))
colormap defult