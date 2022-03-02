close all
clear all
clc

figurePath = 'D:\TuthillLab\figure\220302';

pass_condition_name= {'roi' , 'ball','driver'};
pass_condition_value = {'medial','1','R21D12'};
passCondition(1,1:size(pass_condition_name,2)) = pass_condition_name;
passCondition(2,1:size(pass_condition_value ,2)) = pass_condition_value;

analyse = arcade(figurePath);
[data, metaVariableNames,epochList ]= analyse.eopchListGen(false);
targetID = analyse.dataSelector(epochList,passCondition);
[missFileName,checkedID] = analyse.epochCheck(epochList,targetID);
