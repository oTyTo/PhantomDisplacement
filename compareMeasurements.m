%% compare displacements and calculate average strain
close all;
clear all;
clc;


%% compare displacements
%read data
Path = dir('**/inv_displacement.mat');
for i = 1:length(Path)
   load([Path(i).folder '\' Path(i).name]);
   Y{i} = inv_lastinte;
   X{i} = section;
   sig_length(i) = length(section);
end
for i = 1:length(Path)
    Y{i} = Y{i}(1:min(sig_length));
    X{i} = X{i}(1:min(sig_length));
end
%calcualte area difference
area_first = (sum(Y{1})-0.5*Y{1}(1)-0.5*Y{1}(end))*(X{1}(2)-X{1}(1))*0.43;
for i = 2:length(Path)
    area_diff = (sum(abs(Y{i}-Y{1}))-0.5*abs(Y{i}(1)-Y{1}(1))-0.5*abs(Y{i}(end)-Y{1}(end)))*(X{1}(2)-X{1}(1))*0.43;
    ratio(i-1) = area_diff/abs(area_first);
end


%plot figure
figure
for i = 1:length(Path)
    plot(X{i},Y{i}*0.43);
    if i == 1
        hold on
    end
end
% plot(section,a2,'m')
% plot(section,a3,'r')
% plot(section,a4,'g')
% plot(section,a5,'b')
% plot(loc,dis*1000,'r')
xlabel('signal position (mm)')
ylabel('displacement (um)')
hold off
% title(['Area difference ratio= ' num2str(ratio)])
legend({'1st Measured Result',['2nd Measured Result, similarity = ' num2str(1-ratio(1))],['3rd Measured Result, similarity = ' num2str(1-ratio(2))],['4th Measured Result, similarity = ' num2str(1-ratio(3))],['5nd Measured Result, similarity = ' num2str(1-ratio(4))]},'Location','southwest')
saveas(gcf,'CompareDisplacements.png')

%% calculate strain
for i = 1:length(Path)
    smooth_lastint = smooth(Y{i},0.1,'rloess');
    diff = [];
    for j = 1:1:length(Y{i})-5
        diff(j) = (Y{i}(j+5)-Y{i}(j))*0.43e-3/(X{i}(j+5)-X{i}(j));
    end
    strain(i) = mean(diff);
end
save('avg_strain.mat','strain')
%
% A=[a1;a2;a3;a4;a5];
% f_simu = a5;
% ratio=[];
% for i = 1:4
%     inv_lastinte = A(i,:);
%     area_diff = (sum(abs(inv_lastinte-f_simu))-0.5*abs(inv_lastinte(1)-f_simu(1))-0.5*abs(inv_lastinte(end)-f_simu(end)))*Window_unit;
%     area_sim = (sum(f_simu)-0.5*f_simu(1)-0.5*f_simu(end))*Window_unit;
%     ratio(i) = area_diff/abs(area_sim)
% end
% figure, plot(section,inv_lastinte*0.43,'b'), xlabel('signal position (mm)'),ylabel('displacement (um)');
% hold on

% figure, plot(section(1:end-5),diff,'b'), xlabel('signal position (mm)'),ylabel('strain');
% hold on
% plot(loc,strain,'r')
% xlim([0,3.5])
% hold off
% saveas(gcf,'CompareWithSimulationstrain.png')