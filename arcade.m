classdef arcade
    %   Batch process system for regression analysis
    %   There are four levels of precessing
    %   1. Epoch list: get trail name and meta data
    %   2. Data selector: get the trail id that fits the input condition
    %   3. Data spilt into epoch tree
    %   4. Indivisual computaion: computation on each leaf
    
    properties
        metadataPath ='D:\TuthillLab\video';
        videoPath = 'D:\TuthillLab\video';
        fiugrePath 
        intermediatePath
        
        
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
        
        function [missFileName,checkedID] = epochCheck(obj,epochList,targetID)
            %check and adjust how the eopch change the target.ID list so
            %that it won't reference the video which is not there
            cd(obj.videoPath)
            
            targetList = targetID;
            %epochList = obj.targetList;
            filePassed = [];
            fileFailed = [];
            for i = 1:length(targetList)
                temp_name = epochList(targetList(i)).trail_id;
                fn = [temp_name '_legCam_1.avi'];
                if isfile(fn)
                    
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
        end%end of epochCheck
        
        function epochTree(obj,dataSplit)
            %this function is for data spilting 
        end
        
        
        
        
        
    end%end for methods
end%end for the class
