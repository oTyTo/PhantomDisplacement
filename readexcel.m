% read from excel files
% 
close all;
clear all;
clc;

files = dir('*.csv');
count = zeros(1);
PATH = {'start'};
data = {};
j = 0;
for i = 1:1:length(files)
    NEWPATH = files(i).folder;
    NAME = files(i).name;
    A = xlsread([NEWPATH '/' NAME]);
    if strcmp(NEWPATH , PATH{end})
        count(j) = count(j)+1;
        data{j}{count(j)} = A;        
    else
        j = j+1;
        PATH{j} = NEWPATH;
        count(j) = 1;
        data{j}{count(j)} = A; 
    end
end

%% plot data
loc = data{1}{1}(:,16);
dis = data{1}{1}(:,2);
strain = data{1}{1}(:,8);
figure,
plot(data{1}{1}(:,16),data{1}{1}(:,2))
save('simuwithoutcali.mat','loc','dis','strain')