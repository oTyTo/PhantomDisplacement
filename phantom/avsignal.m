%% average of signals 
% 
close all;
clear all;
clc;
tic


%% some parameters
num_f = 3;      %num of files per measurement

%% load data from .mat
files = dir('**/*.mat'); 

count = zeros(1);
PATH = {'start'};
data = {};
startTime = [];
intevalTime = [];
j = 0;
for i = 1:1:length(files)
    NEWPATH = files(i).folder;
    NAME = files(i).name;
    if strcmp(NAME,'simudis.mat')==0
        load([NEWPATH '/' NAME]);
        if strcmp(NEWPATH , PATH{end})
            count(j) = count(j)+1;
            data{j}(count(j),:) = A;
            startTime(j,count(j)) = Tstart;
            intevalTime(j,count(j)) = Tinterval;
        else
            %check if previous data directory have 32 groups of signals
            if j>0
                if count(j)<32
                    disp([num2str(j) 'th group has only ' num2str(count(j)) ' of signals'])
                end
            end
            j = j+1;
            PATH{j} = NEWPATH;
            count(j) = 1;
            data{j}(count(j),:) = A;
            startTime(j,count(j)) = Tstart;
            intevalTime(j,count(j)) = Tinterval;
        end
    end
end
toc

%% check time inteval and align data based on time
if length(unique(intevalTime))>1
    disp([num2str(length(unique(intevalTime))) ' signals have different time inteval'])
end
if length(unique(startTime))>1
    disp([num2str(length(startTime(startTime ~= median(median(startTime))))) ' signals have different start time'])
end
if sum(startTime>0)>0
    disp('Some signals have positive start time')
end
%data alignment, assume all start times are negative
cali_data = {};
[M,N] = size(startTime);
for i = 1:1:M
    for j = 1:1:count(i)
        point0 = ceil(abs(startTime(i,j))/intevalTime(i,j));
        %point0 = point0 - 1; %when use pico 6000, mat file has one less data point
         tempdata = [data{i}(j,point0+1:end),zeros(1,point0)];
         cali_data{i}(j,:) = tempdata(1:25000);               %set all signal to same length
    end
end

%% average of each sample group(num_f folders)
avg = [];
for i = 1:num_f:length(count)
    tempsum = zeros(1,size(cali_data{i},2));
    for j = 0:1:num_f-1
        tempsum = tempsum + median(cali_data{i+j});
    end
    avg(floor((i-1)/num_f)+1,:) = tempsum./num_f;
end
save('avdata.mat','avg');
toc

%% load data from .csv files
% files = dir('**/*.csv');
% count = zeros(1);
% PATH = {'start'};
% data = {};
% j = 0;
% for i = 1:1:length(files)
%     NEWPATH = files(i).folder;
%     NAME = files(i).name;
%     A = xlsread([NEWPATH '/' NAME]);
%     if strcmp(NEWPATH , PATH{end})
%         count(j) = count(j)+1;
%         data{j}{count(j)} = A;        
%     else
%         j = j+1;
%         PATH{j} = NEWPATH;
%         count(j) = 1;
%         data{j}{count(j)} = A; 
%     end
% end
% toc
% %% calibrate so that every sequence start from 0 time
% cali_data = {};
% for i = 1:1:length(data)
%     for j = 1:1:length(data{i})
%         [rows,~] = find(data{i}{j}(:,1)<0);
%         cali_data{i}(j,:) = [data{i}{j}(length(rows)+1:end,2);zeros(length(rows),1)];
%     end
% end
% 
% %% average of each sample group(3 folders)
% num_f = 3;
% avg = [];
% for i = 1:num_f:length(count)
%     tempsum = zeros(1,size(cali_data{i},2));
%     for j = 0:1:num_f-1
%         tempsum = tempsum + mean(cali_data{i+j});
%     end
%     avg(floor(i/num_f)+1,:) = tempsum./num_f;
% end
% save('avdata.mat','avg');


