classdef arcade
    %   Batch process system for regression analysis
    %   There are four levels of precessing
    %   1. Epoch list: get trail name and meta data
    %   2. Data selector: get the trail id that fits the input condition
    %   3. Data spilt into epoch tree
    %   4. Indivisual computaion: computation on each leaf
    
    properties
        metadataPath ='D:\TuthillLab\video';
        videoPath = 'C:\Users\USER\Downloads\Data';
        fiugrePath
        intermediatePath
        
        
        
    end
    
    properties (Access = private)
        data%data table file
        metaVariableNames%metadate name
        epochList%list of unique epochs defined by trail name
        targetID%ID of the trails after the selector
        missFileName%files you need to down load
        checkedID%it will only reference the video on the computer
    end
    
    
    methods
        function obj = arcade(fiugrePath)
            if nargin ==1
                obj.fiugrePath = fiugrePath;
            end
        end
        
        %there are bugs on setting attribute, fix later;
        
        
        function [data, metaVariableNames,epochList ]= eopchListGen(obj,makeNewList)
            %For generating epoch list for data selector, save and
            %re-organize the data
            
            %input obj.metadataPath, location of data.csv or data.mat
            %input makeNewList = ture, rerun the list when you add new
            %datas; false, load old processed data.
            
            %output data metaVaribleNames
            %assign the data and metaVariableNames properity
            if makeNewList == true %make new list and load data file
                %load
                cd(obj.metadataPath)
                % Initialize variables.
                filename = 'D:\TuthillLab\data\data.csv';
                delimiter = ',';
                startRow = 2;
                formatSpec = '%C%C%C%C%C%q%C%q%f%f%f%C%f%f%f%f%f%[^\n\r]';
                
                % Open the text file.
                fileID = fopen(filename,'r');
                
                % Read columns of data according to the format.
                % This call is based on the structure of the file used to generate this
                % code. If an error occurs for a different file, try regenerating the code
                % from the Import Tool.
                dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                % Close the text file.
                fclose(fileID);
                % Post processing for unimportable data.
                % No unimportable data rules were applied during the import, so no post
                % processing code is included. To generate code which works for
                % unimportable data, select unimportable cells in a file and regenerate the
                % script
                % Create output variable: save and assign to obj
                %give direct output, easier to debug
                data = table(dataArray{1:end-1}, 'VariableNames', {'roi','driver','ball','animal','trial','frame','activity','prediction','L1C_flex','L2C_flex','L3C_flex','analyze','L1_rest','L1_walk','L1_groom','L1_other','L3_groom'});
                metaVariableNames = {'roi','driver','ball','animal','trial','frame','activity','prediction','L1C_flex','L2C_flex','L3C_flex','analyze','L1_rest','L1_walk','L1_groom','L1_other','L3_groom'};
                obj.data = data;
                obj.metaVariableNames = metaVariableNames;
                
                
                
                save('data','data') %save and overwrite the old one
                save('metaVariableNames','metaVariableNames')
                % Clear temporary variables
                clearvars filename delimiter startRow formatSpec fileID dataArray ans;
                %%generate epoch list
                trialz = data(:,5);
                trialu = unique(trialz);%unique trail id
                %extract the meta data;
                
                parfor trail_idx = 1:size(trialu ,1)
                    epochList(trail_idx).trail_id = char(trialu{trail_idx,1});
                    %extract the experiment parameter
                    for i = 1:3000:size(data,1)%loop over whole data file;
                        if trialz{i,1}== trialu{trail_idx,1}
                            epochList(trail_idx).roi =char(data{i,1});
                            epochList(trail_idx).driver = char(data{i,2});
                            epochList(trail_idx).ball = char(data{i,3});
                            epochList(trail_idx).animal= char((data{i,4}));
                            break
                        end
                    end
                end
                
                %add a warining if there is a skip of condition
                if length(epochList)~=size(trialu,1)
                    disp('some trails are missing, change step size')
                end
                %add other traces later
                obj.epochList = epochList;
                save('epochList','epochList')
                
            elseif makeNewList == false
                cd(obj.metadataPath)
                load('data.mat')
                load('metaVariableNames')
                load('epochList')
                obj.data = data;
                obj.metaVariableNames = metaVariableNames;
                obj.epochList = epochList;
            end
        end %end of epochListGen
        
        function targetID = dataSelector(obj,epochList,passCondition)
            
            %input: the condition we want to select
            %output: the trail ID for such trails
            %tempList = obj.epochList;
            tempList = epochList;
            traget_trail = [];%should be target
            for i = 1:size(passCondition,2)
                if isempty(traget_trail )
                    for j = 1:length(tempList)%index of each trail
                        temp_match = passCondition{2,i};
                        if strcmp(tempList(j).(passCondition{1,i}),passCondition{2,i})
                            traget_trail = [ traget_trail j];
                        end
                    end
                else
                    ph = [];
                    for k = 1:length(traget_trail)
                        s= traget_trail(k);
                        if strcmp(tempList(s).(passCondition{1,i}),passCondition{2,i})
                            ph = [ph s];
                        end
                    end
                    traget_trail = ph;
                    
                end
            end
            targetID = traget_trail;
        end %end of data selector
        
        function [missFileName,checkedID,trailList] = epochCheck(obj,epochList,targetID)
            %check and adjust how the eopch change the target.ID list so
            %that it won't reference the video which is not there
            cd(obj.videoPath)
            
            targetList = targetID;
            %epochList = obj.targetList;
            filePassed = [];
            fileFailed = [];
            trailList = [];
            for i = 1:length(targetList)
                temp_name = epochList(targetList(i)).trail_id;
                fn = [temp_name '_legCam_1.avi'];
                if isfile(fn)
                    trailList = [ trailList; temp_name];
                    filePassed = [filePassed; targetList(i)];
                else
                    
                    fileFailed = [fileFailed; targetList(i)];
                end
            end
            
            
            
            if ~isempty(fileFailed)
                missFileName =[];
                for j = 1:length(fileFailed)
                    missFileName = [missFileName; epochList(fileFailed(j)).trail_id];
                end
                missFileName = missFileName;
            end
            
            missFileName = [];
            checkedID = filePassed;
            cd(obj.metadataPath)
            save('trailList','trailList')
        end%end of epochCheck
        
        function epochTree(obj,dataSplit)
            %this function is for data spilting
        end
        
        function [gluedDataMatrix,gluedMetaData,douball_mask,ref_frame_m,down_test_frame]= videoReaderGlue(obj,trailList,data)
            tic
            sub_fn = '_legCam_1.avi';
            video_temp = [];
            calcuim_trace_temp = [];
            prediction_temp = [];
            L1C_f_temp = [];
            L2C_f_temp = [];
            L3C_f_temp = [];
            L1_rest_temp = [];
            L1_walk_temp = [];
            L1_groom_temp = [];
            L1_other_temp = [];
            L3_groom_temp = [];
            
            
            for tid = 1:size(trailList,1)%differnt epoch from same animal
                temp_name = trailList(tid,:);
                nov = [temp_name sub_fn]; %file name
                trail_id =nov(1:18); % trail ID in the data table
                %extract the meta data for that video
                trail_list = string(data{:,5});% list of trails extracted from the data
                trail_target_index = find(trail_list == trail_id);% index of the video we want to analyze
                sub_data = data(trail_target_index,:); %extract the subset of the data that is realted to the video trail
                
                
                frame_number_pre = double(string(sub_data{:,6}));%frame number for alignment
                calcuim_trace_pre = double(string(sub_data{:,7}));%calcuim traces
                prediction_pre = double(string(sub_data{:,8}));
                L1C_f_pre = double(string(sub_data{:,9}));
                L2C_f_pre = double(string(sub_data{:,10}));
                L3C_f_pre = double(string(sub_data{:,11}));
                
                analyze_chek = double(string(sub_data{:,12}));%if zero then don't analyse
                
                L1_rest_pre = double(string(sub_data{:,13}));
                L1_walk_pre = double(string(sub_data{:,14}));
                L1_groom_pre = double(string(sub_data{:,15}));
                L1_other_pre = double(string(sub_data{:,16}));
                L3_groom_pre = double(string(sub_data{:,17}));
                
                flag_pole = analyze_chek;%.*(walk_flag==0);%this is and opertaion
                
                calcuim_trace = calcuim_trace_pre(flag_pole==1);%clean the calcuim trace
                prediction_trace  = prediction_pre(flag_pole==1);%clean up the prediction
                frame_number = frame_number_pre(flag_pole==1);%clean the frame number
                L1C_f = L1C_f_pre(flag_pole==1);
                L2C_f = L2C_f_pre(flag_pole==1);
                L3C_f = L3C_f_pre(flag_pole==1);
                L1_rest = L1_rest_pre(flag_pole==1);
                L1_walk = L1_walk_pre(flag_pole==1);
                L1_groom = L1_groom_pre(flag_pole==1);
                L1_other = L1_other_pre(flag_pole==1);
                L3_groom = L3_groom_pre(flag_pole==1);
                
                %video read
                cd(obj.videoPath)
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
                
                %gluing starts here
                video_temp = [video_temp, video_matrix];
                calcuim_trace_temp = [calcuim_trace_temp; calcuim_trace];
                prediction_temp = [prediction_temp; prediction_trace];
                L1C_f_temp = [ L1C_f_temp;L1C_f];
                L2C_f_temp = [L2C_f_temp;L2C_f ];
                L3C_f_temp = [L3C_f_temp;L3C_f];
                L1_rest_temp = [L1_rest_temp;L1_rest];
                L1_walk_temp = [L1_walk_temp;L1_walk ];
                L1_groom_temp = [L1_groom_temp;L1_groom];
                L1_other_temp = [L1_other_temp;L1_other];
                L3_groom_temp = [L3_groom_temp;L3_groom];
            end
            gluedDataMatrix = video_temp ;
            gluedMetaData.calcuim_trace = calcuim_trace_temp;
            gluedMetaData.prediction = prediction_temp;
            gluedMetaData.L1C_f= L1C_f_temp;
            gluedMetaData.L2C_f = L2C_f_temp;
            gluedMetaData.L3C_f = L3C_f_temp;
            gluedMetaData.L1_rest =  L1_rest_temp;
            gluedMetaData.L1_walk = L1_walk_temp;
            gluedMetaData.L1_groom = L1_groom_temp;
            gluedMetaData.L1_other = L1_other_temp;
            gluedMetaData.L3_groom = L3_groom_temp;
            toc
        end%end of video reader glue
        
        function [V_r,S_r,A,U] = videoPC(obj,k_components,dataMatrix)
            tic
            [U,S,V] = svds(dataMatrix,k_components);%SVD of data matrix
            toc
            sv = diag(S);
            figure
            stairs((cumsum(sv)/sum(sv)).*(cumsum(sv)/sum(sv)))
            title('Explained varience')
            V_conj = V';
            r =k_components;%take first r pc
            V_r = V_conj(1:r,:);%reduce to r-dimensions
            S_r = S(1:r,1:r);%reduce to r-dimensions
            A = S_r*V_r;%project to r-dimensions
           
        end%end of video PC
        
        function[B,FitInfo,imr,y_hat]= lassoReconstruct(obj, U,A,targertTrace,r)
             % lasso
            [B,FitInfo] = lasso(A',targertTrace,'Alpha',0.5);%elsAtic net
            [~,min_idx]= min(FitInfo.MSE);%find the minium lambd
            singular_value_rank = min_idx; %use the correspond singular value
            Bb= B(:,singular_value_rank)';%get n-th row out
            U_r = U(:,1:r);%reduce U to r by b
            imr = U_r*Bb';
            y_hat = (Bb*A)+FitInfo.Intercept(singular_value_rank )*ones(1,length((Bb*A)));
        end%end of lasso reconstruct 
        
        
    end%end for methods
end%end for the class
