%% compare displacements and calculate average strain
close all;
clear all;
clc;


%% compare displacements
%read data
Path = dir('**/inv_displacement.mat');
%%
for i = 1:length(Path)
   load([Path(i).folder '\' Path(i).name]);
   Y{i} = inv_lastinte;
   X{i} = section;
   sig_length(i) = length(section);
end
%%
for i = 1:length(Path)
    Y{i} = Y{i}(1:min(sig_length));
    X{i} = X{i}(1:min(sig_length));
end
%avg in groups
median_dis(1,:) = median([Y{1};Y{2};Y{3};Y{4};Y{5}]);
median_dis(2,:) = median([Y{6};Y{7};Y{8};Y{9}]);
median_dis(3,:) = median([Y{10};Y{11};Y{12};Y{13};Y{14}]);
median_dis(4,:) = median([Y{15};Y{16};Y{17};Y{18};Y{19}]);

%plot figure
figure
% for i = 1:4
%     plot(X{i},median_dis(i,:)*0.43);
%     if i == 1
%         hold on
%     end
% end
plot(X{1},median_dis(4,:),'m')
hold on
plot(X{1},median_dis(2,:),'r')
plot(X{1},median_dis(1,:),'g')
plot(X{1},median_dis(3,:),'b')
% plot(loc,dis*1000,'r')
xlabel('signal position (mm)')
ylabel('displacement (um)')
hold off
% title(['Area difference ratio= ' num2str(ratio)])
legend({'10.3kPa','29.0kPa','31.1kPa','48.6kPa'},'Location','southwest')
saveas(gcf,'CompareDisplacements.png')

% %% calculate strain
% for i = 1:length(Path)
%     smooth_lastint = smooth(Y{i},0.1,'rloess');
%     diff = [];
%     for j = 1:1:length(Y{i})-5
%         diff(j) = (Y{i}(j+5)-Y{i}(j))*0.43e-3/(X{i}(j+5)-X{i}(j));
%     end
%     strain(i) = mean(diff);
% end
% save('avg_strain.mat','strain')
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