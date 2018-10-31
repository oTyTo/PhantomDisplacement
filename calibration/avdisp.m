%% process using average signal
% with pure cover
close all;
clear all;
clc;

mkdir signal 
mkdir displacements
load('avdata.mat')
% for i = 1:1:size(avg,1)
%     figure
%     plot (avg(i,:));
% end
% %dont use pure cover data
% avg = avg(1:end-1,:);

%% down sample the input siganl
% avg= downsample(avg',4)';
%% deal with infi and nan 
%set infi, nan to max
avg(isinf(avg)) = 1;
avg(isnan(avg)) = 1;
%set nan to []
% [M, N] = find(isnan(avg));
% nanrows = unique(M);
% for i = 1:1:length(nanrows)
%     avg(nanrows(i)-(i-1),:)=[];
% end
%display first and last signal to choose a proper parameter
figure
plot(avg(1,:)); hold on;
plot(avg(1,:)); hold off;
legend('ref','com');    
%4671-17950
%% fft of original signal
% sig = avg(1,:);
% L = length(sig);
% t = 2e-5;
% fs = L/t;
% Y = fftshift(fft(sig));
% SP = abs(Y/L);
% f = 1/t*(-L/2:L/2-1);
% 
% %apply a bandpass filter
% [z,p,k] = butter(10,[4e6 9e6]*2/fs,'bandpass');
% sos = zp2sos(z,p,k);
% bandsig = sosfilt(sos,sig);
% Y = fftshift(fft(bandsig));
% bandSP = abs(Y/L);
% f = 1/t*(-L/2:L/2-1);
% figure
% subplot(2,1,1)
% plot (f,SP)
% subplot(2,1,2)
% plot(f,bandSP);
% figure
% subplot(2,1,1)
% plot (sig)
% subplot(2,1,2)
% plot(bandsig);
% %% bandpass of original signal
% L = size(avg,2);
% t = 2e-5;
% fs = L/t;
% %apply a bandpass filter
% [z,p,k] = butter(10,[4e6 9e6]*2/fs,'bandpass');
% sos = zp2sos(z,p,k);
% avg = sosfilt(sos,avg,2);



%% set all parameters;
%0.6mm,8mm, 2.1/2.2
%start:3292 
%first reflect: 4754-5054
method        = 2;                       % 1 is with high pass filter, 2 no filter
para.window   = 30;                     % window size
para.delt_w   = round(para.window/5);   % window overlap
para.tau      = 50;                      % search range
para.startP   = 4356;                    % selected starting point
para.endP     = 4596;                   % selected ending point
para.fs       = 1.25e9;                   % sampling rate of your oscilloscope
para.cut_freq = 2.5e6;                   % remove low frequency motion
para.order    = 4;                       % order of the filter

%% motion estimation
disp_matrix = [];
disp_accum = [];
corr_matrix = [];
for i = 1:(size(avg,1)-1)
    ref = avg(i,:);
    com = avg(i+1,:);
    [filt_ref, filt_com, displacement] = motionEst(ref,com,para,method);
    f = figure();
    plot(filt_ref);
    hold on,plot(filt_com,'r');

    saveas(f,['signal\' num2str((i-1)*5) '-' num2str(i*5) '.png'])
%     savefig(f,['signal\' num2str((i-1)*5) '-' num2str(i*5) '.fig'])
    close(f);
    disp_matrix = [disp_matrix; displacement(1,:)];
    disp_accum = [disp_accum; sum(disp_matrix,1)];
    corr_matrix = [corr_matrix; displacement(2,:)];

    if isempty(displacement) == 1
        disp('your input parameters are wrong');
    else
        frame = size(displacement,2);            % 1st row is displacement idx, 2st row is it correlation


        % you can convert your moving point to real distance
        dist_point_convert = 10/7000;  % mm per each point
        x_axis = linspace(7.88,2.81,frame);
        disp_curve = displacement * dist_point_convert;
    end
end

%% visualize original displacements
for i = 1:1:size(disp_matrix,1)

    g = figure();
    subplot(2,1,1);
    plot(disp_matrix(i,:));       % displacement curve, unit: point
    subplot(2,1,2);
    plot(corr_matrix(i,:));       % correlation accuracy
    saveas(g,['displacements\' num2str((i-1)*5) '-' num2str(i*5) '.png'])
    %savefig(f,['signal\' num2str((i-1)*5) '-' num2str(i*5) '.fig'])
    close(g);
end

%% process results
%use rows with mean correlation >0.9
% rows = find(mean(corr_matrix,2)>0.9);
% tempm = disp_matrix(rows,:);
tempm = disp_matrix;
% change all value>0 or abs>50 to mean
[m,n] = find(tempm>0 | tempm<-5);
for i = 1:1:length(m)
    tempm(m(i),n(i)) = median(tempm(tempm(:,n(i))<0 & tempm(:,n(i))>-5,n(i)));
end
tempacc = zeros(size(tempm));
for i = 1:1:size(tempm,1)
    tempacc(i,:) = sum(tempm(1:i,:),1);
end
avg_dis  = mean(tempacc(end,:),2);
avg_diff = avg_dis*0.43e-3/0.6/1.4;

%% get mean displacement accross the time interval
%mean_dis = mean(disp_accum(end,:),2);
[filt_ref, filt_com, displacement] = motionEst(avg(1,:),avg(end,:),para,method);
t = para.startP:1:para.endP;
t = t*8e-10;
figure();
plot(t,filt_ref);
hold on,plot(t,filt_com,'r');
title('First and last signal'),xlabel('time (s)'),ylabel('voltage (V)')
figure();
plot(displacement(1,:)),title('displacement'),ylim([-100,0]);       % displacement curve, unit: point
avg_dis = mean(displacement,2);
avg_diff = avg_dis*0.43e-3/0.6/1.4;
saveas(gcf,['upperstrain=' num2str(avg_diff(1)) '.png'])

%% visualize results
section = 0:(9.32/16300*10400/size(displacement,2)):9.32/16300*10400-(9.32/16300*10400/size(displacement,2));
%integrated displacement
%figure, plot(section,tempm(end,:)*0.43), xlabel('signal position (mm)'),ylabel('displacement (um)');
inv_lastinte = tempacc(end,end)-flip(tempacc(end,:));
figure, plot(section,inv_lastinte*0.43,'.b'), xlabel('signal position (mm)'),ylabel('displacement (um)');
p = polyfit(section,inv_lastinte*0.43,3);
f = polyval(p,section);
hold on
plot(section,f,'--r')
saveas(gcf,['integrated displacement' num2str(0) '-' num2str(19) '.png'])
%real displacement
[filt_ref, filt_com, displacement] = motionEst(avg(1,:),avg(end,:),para,method);
inv_last = displacement(1,end)-flip(displacement(1,:));
figure, plot(section, inv_last*0.43),xlabel('signal position (mm)'),ylabel('displacement (um)');
saveas(gcf,['cross displacement' num2str(0) '-' num2str(19) '.png'])
%2d integral
figure
imagesc(tempacc),colorbar
xticklabels({'0.838','1.676','2.514','3.352','4.190','5.028','5.866'})
xlabel('signal position (mm)')
saveas(gcf,'Integraded2D.png')
%2d difference to 1
figure
imagesc(tempm),colorbar
xticklabels({'0.838','1.676','2.514','3.352','4.190','5.028','5.866'})
xlabel('signal position (mm)')
saveas(gcf,'NeighborDisplacements2D.png')

% %%
% rows = find(mean(corr_matrix,2)>0.9);
% for i = 1:1:length(rows)
%     figure
%     plot(avg(rows(i),:)); hold on;
%     plot(avg(rows(i)+1,:));hold off;
% end