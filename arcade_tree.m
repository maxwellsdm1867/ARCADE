%acarade epoch tree
%use to pull batch of trail name out for a given conditon 
close all 
clear all
clc
%% Specify working directory and load the raw and metadata
cd('D:\TuthillLab\data')%go to the data folder, it should contain chris_data.mat and the video you wish to analyze
load('D:\TuthillLab\data\chris_data.mat')%loading the convterd csv table, this contains the metadata and calcuim trace
data = data(2:end,:);%take out the first row, since origianlly it's just text
%%get the trail name for the specific condition 
roiz = data(:,1);
driverz = data(:,2);
ballz = data(:,3);
animalz = data(:,4);
trialz = data(:,5);
%get the unique roi, dirver, animal, ball and animal 
roiu = unique(roiz);
driveru = unique(driverz);
ballu = unique(ballz);
animalu = unique(animalz);
tic
traillist = [];
parfor i = 1:size(data,1)
    
    if driverz{i,1}=='R21D12' && ballz{i,1} == '1' && roiz{i,1}== roiu{1,1}
         traillist = [ traillist; trialz{i,1} ];
    end
end
toc

traillist = unique(traillist );

%to be continued