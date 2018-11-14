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

%% compare with simulations
load('simu.mat')
%calculate lines
p_simu = constrainfit(loc(1:1000)',dis(1:1000)'*1000,0,0,10);
f_simu = polyval(p_simu,X{1});
% g = corr2(inv_lastinte*0.43,f_simu);
% p_line = constrainfit(inv_lastinte*0.43,f_simu,0,0,1);
% line = polyval(p_line,inv_lastinte*0.43);
%calcualte area difference
area_simu = (sum(f_simu)-0.5*f_simu(1)-0.5*f_simu(end))*(X{1}(2)-X{1}(1));
for i = 1:length(Path)
    area_diff = (sum(abs(Y{i}*0.43-f_simu))-0.5*abs(Y{i}(1)*0.43-f_simu(1))-0.5*abs(Y{i}(end)*0.43-f_simu(end)))*(X{1}(2)-X{1}(1));
    ratio(i) = area_diff/abs(area_simu);
end
%plot
figure
plot(X{1},f_simu)
hold on
for i = 1:length(Path)
    plot(X{i},Y{i}*0.43);
end
hold off
xlabel('signal position (mm)')
ylabel('displacement (um)')
legend({'Simulation Result',['1st Measured Result, similarity = ' num2str(1-ratio(1))],['2nd Measured Result, similarity = ' num2str(1-ratio(2))],['3rd Measured Result, similarity = ' num2str(1-ratio(3))],['4th Measured Result, similarity = ' num2str(1-ratio(4))],['5nd Measured Result, similarity = ' num2str(1-ratio(5))]},'Location','southwest')
saveas(gcf,'CompareWithSimulation.png')

%% compare mean and simu
median_data = median([Y{1};Y{2};Y{3};Y{4};Y{5}]);
%calculate difference
area_diff = (sum(abs(median_data*0.43-f_simu))-0.5*abs(median_data(1)*0.43-f_simu(1))-0.5*abs(median_data(end)*0.43-f_simu(end)))*(X{1}(2)-X{1}(1));
ratio = area_diff/abs(area_simu);
figure
plot(X{1},f_simu)
hold on
plot(X{1},median_data*0.43)
hold off
xlabel('signal position (mm)')
ylabel('displacement (um)')
legend({'Simulation Result',['Avegare of measured data, similarity = ' num2str(1-ratio)]},'Location','southwest')
saveas(gcf,'MeanWithSimulation.png')

% figure, plot(inv_lastinte*0.43,f_simu),xlabel('experiment displacement (mm)'),ylabel('simulation displacement (um)');
% hold on
% plot(inv_lastinte*0.43,line,'--')
% title(['Correlation Coefficient= ' num2str(g)])
% hold off
% saveas(gcf,'CC.png')
% 
% figure, plot(section,inv_lastinte*0.43,'b'), xlabel('signal position (mm)'),ylabel('displacement (um)');
% hold on
% plot(loc,dis*1000,'r')
% xlim([0,3.5])
% hold off
% title(['Area difference ratio= ' num2str(ratio)])
% legend({'Measured Result','Simulation Result'},'Location','southwest')
% saveas(gcf,'CompareWithSimulation.png')
% 
% figure, plot(section(1:end-5),diff,'b'), xlabel('signal position (mm)'),ylabel('strain');
% hold on
% plot(loc,strain,'r')
% xlim([0,3.5])
% hold off
% saveas(gcf,'CompareWithSimulationstrain.png')