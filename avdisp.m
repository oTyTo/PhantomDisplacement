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
%% phantom specs
P_height = 9.02;%(mm)
P_cali = 4356;%(points)
P_bottom = 18060;
P_unit = P_height/(P_bottom-P_cali);

%% down sample the input siganl
% avg= downsample(avg',4)';
%% deal with infi and nan 
%set infi, nan to 2, 1
avg(isinf(avg)) = 0.1;
avg(isnan(avg)) = 0.1;
%set nan to []
% [M, N] = find(isnan(avg));
% nanrows = unique(M);
% for i = 1:1:length(nanrows)
%     avg(nanrows(i)-(i-1),:)=[];
% end
%display first and last signal to choose a proper parameter
figure
subplot(2,1,1)
plot(avg(1,:)); 
subplot(2,1,2)
plot(avg(end,:),'r');


%% more phantom specs (points)
P_top = 8292;
P_top_press = 7749;
P_bot = 17950;
P_bot_press = 17320;
%% fft of original signal
%
sig = avg(1:end,P_top:P_bot_press);
dt_time = 8e-10;
freq_sam = 1/dt_time;
NFFT = size(sig,2);% depending on the length of the time signal
freq_reso = freq_sam/NFFT;
freq = 0:freq_reso:freq_sam/2-freq_reso;
freq = freq/1e6;

%
wav_filt_fft = abs(fft(sig,NFFT,2)./NFFT);
% wav_filt_fft = smooth(wav_filt_fft,20);
figure
plot(freq,(wav_filt_fft(1,1:NFFT/2)),'r')
xlim([0 15])
bandsig = sig;
%% apply band pass filter 
%apply a bandpass filter
[z,p,k] = butter(2,[7e6 11e6]*2/freq_sam,'bandpass');
sos = zp2sos(z,p,k);
bandsig = sosfilt(sos,sig,2);
Y = abs(fft(bandsig,NFFT,2)./NFFT);
figure
subplot(2,1,1)
plot(freq,(wav_filt_fft(1,1:NFFT/2)),'r')
xlim([0 15])
subplot(2,1,2)
plot(freq,(Y(1,1:NFFT/2)),'r')
xlim([0 15])
figure
subplot(2,1,1)
plot (sig(1,:))
subplot(2,1,2)
plot(bandsig(1,:));




%% set all parameters;
%0.6mm,8mm, 2.1/2.2
%start:2529 
%first reflect: 3659-3890
%end : 15310
method        = 2;                       % 1 is with high pass filter, 2 no filter
para.window   = 1500;                     % window size
para.delt_w   = round(para.window/50);   % window overlap
para.tau      = 100;                      % search range
para.startP   = 1;                    % selected starting point
para.endP     = size(bandsig,2);                   % selected ending point
para.fs       = 1.25e9;                   % sampling rate of your oscilloscope
para.cut_freq = 2.5e6;                   % remove low frequency motion
para.order    = 4;                       % order of the filter

%% motion estimation
disp_matrix = [];
disp_accum = [];
corr_matrix = [];
for i = 1:(size(bandsig,1)-1)
    ref = bandsig(i,:);
    com = bandsig(i+1,:);
    [filt_ref, filt_com, displacement] = motionEst(ref,com,para,method);
    f = figure();
    subplot(2,1,1)
    plot(filt_ref);
    hold on,plot(filt_com,'r'),ylim([-0.02,0.02]);
    subplot(2,1,2)
    plot(displacement(1,:));

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
tempm = disp_matrix;
%change data points with corrlation <0.8 to mean
[m,n] = find(corr_matrix<0.7);
for i = 1:1:length(m)
    tempm(m(i),n(i)) = mean(tempm(corr_matrix(:,n(i))>=0.7,n(i)));
end

% change all value>0 or abs>30 to mean
[m,n] = find(tempm>0 | tempm<-40);
for i = 1:1:length(m)
    tempm(m(i),n(i)) = mean(tempm(tempm(:,n(i))<0 & tempm(:,n(i))>-40,n(i)));
end
tempacc = zeros(size(tempm));
for i = 1:1:size(tempm,1)
    tempacc(i,:) = sum(tempm(1:i,:),1);
end


%% visualize results
P_select_length = P_unit*(para.endP-para.startP);
Window_unit = P_select_length/size(displacement,2);
section = 0:Window_unit:P_select_length-Window_unit;
%integrated displacement
%figure, plot(section,tempm(end,:)*0.43), xlabel('signal position (mm)'),ylabel('displacement (um)');
inv_lastinte = tempacc(end,end)-flip(tempacc(end,:));
figure, plot(section,inv_lastinte*0.43,'.b'), xlabel('signal position (mm)'),ylabel('displacement (um)');
p = constrainfit(section,inv_lastinte*0.43,0,0,10);
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
%% calculate strain
smooth_lastint = smooth(inv_lastinte,0.1,'rloess');
diff = [];
diff_ori = [];
diff_smooth = [];
for i = 1:1:length(f)-5
    diff(i) = (f(i+5)-f(i))*0.43e-3/(5*Window_unit);
    diff_ori(i) = (inv_lastinte(i+5)-inv_lastinte(i))*0.43e-3/(5*Window_unit);
    smooth_diff_ori = smooth(diff_ori,7);
    diff_smooth(i) = (smooth_lastint(i+5)-smooth_lastint(i))*0.43e-3/(5*Window_unit);
end
%diff(diff>0) = median(diff,2);
mean_diff = mean(diff,2);
figure
subplot(2,1,1),plot(section,f),title('displacement'),xlabel('thickness (mm)'),ylabel('displacement (um)')
subplot(2,1,2),plot(section(1:end-5),diff),title('strain'),xlabel('thickness (mm)'),ylabel('strain')
saveas(gcf,['HighOrderFitting' num2str(mean_diff) '.png'])
figure
subplot(2,1,1),plot(section,inv_lastinte),title('displacement'),xlabel('thickness (mm)'),ylabel('displacement (um)')
subplot(2,1,2),plot(section(1:end-5),diff_ori),title('strain'),xlabel('thickness (mm)'),ylabel('strain')
saveas(gcf,['OriginalData' num2str(mean_diff) '.png'])
figure
subplot(2,1,1),plot(section,inv_lastinte),title('displacement'),xlabel('thickness (mm)'),ylabel('displacement (um)')
subplot(2,1,2),plot(section(1:end-5),smooth_diff_ori),title('strain'),xlabel('thickness (mm)'),ylabel('strain')
saveas(gcf,['SmoothedOriginal' num2str(mean_diff) '.png'])
figure
subplot(2,1,1),plot(section,smooth_lastint),title('displacement'),xlabel('thickness (mm)'),ylabel('displacement (um)')
subplot(2,1,2),plot(section(1:end-5),diff_smooth),title('strain'),xlabel('thickness (mm)'),ylabel('strain')
saveas(gcf,['DiscreteFitting' num2str(mean_diff) '.png'])
%%
% rows = find(mean(corr_matrix,2)>0.9);
% for i = 1:1:length(rows)
%     figure
%     plot(avg(rows(i),:)); hold on;
%     plot(avg(rows(i)+1,:));hold off;
% end
% load('simudis.mat')
% figure, plot(section,inv_lastinte*0.43,'b'), xlabel('signal position (mm)'),ylabel('displacement (um)');
% hold on
% plot(loc,dis*1000,'r')
% legend({'Measured Result','Simulation Result'},'Location','southwest')
% saveas(gcf,'CompareWithSimulation.png')