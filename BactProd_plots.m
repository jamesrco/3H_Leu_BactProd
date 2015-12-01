% BLATZ_VICE_BactProd_plots.m
% 11/3/2013 JRC
% plots & workups of bacterial production data from BLATZ II and NA VICE
% cruises (KN207-1 and KN207-3)

% presumes you have already performed all necessary averaging and calculations using
% plotQuenched.m and riBPdata.m

% riBPdata.m should have generated the following files as output:
%   BLATZ_VICE_BactProd_results_by_sample_ID.csv
%   BLATZ_VICE_BactProd_Station_data_by_station_no.csv
%   BLATZ_VICE_BactProd_results_replicate_averaged_by_CTD.csv
%   BLATZ_VICE_BactProd_results_replicate_averaged_by_experiment.csv
%   BLATZ_VICE_BactProd_calcs_20xxxxxx.mat, where xxxxxx is the date

clear all; close all;

scrsz = get(0,'ScreenSize'); % define a screen size variable so I can make figures that look decent

%% define values of some constants

ID = 1; % isotope dilution; value of 1 will yield most conservative estimate
leu_BPAA = 0.073; % mol fraction leucine : total bacterial amino acids, Simon & Azam 1989
leu_BPAA_uncert = 0.0191; % associated uncertainty, Simon & Azam 1989
C_prot = 0.86; % mass ratio cellular carbon : protein for marine bacteria, Simon & Azam 1989
MW_leu = 131.2; % molecular weight of leucine

%% load data

% BP data first

BP_replicate_averaged_by_CTD = csvread('BLATZ_VICE_BactProd_results_replicate_averaged_by_CTD.csv',1,0);

size = size(BP_replicate_averaged_by_CTD);
numcols = size(2);

fid = fopen('BLATZ_VICE_BactProd_results_replicate_averaged_by_CTD.csv');

BP_replicate_averaged_by_CTD_headers = textscan(fid,'%s',numcols,'delimiter',',');

fclose(fid);

% extract data from matrix into separate vectors and then assign vector
% names

for i=1:numcols
    colname = BP_replicate_averaged_by_CTD_headers{1,1}(i);
    colname = char(colname);
    colnumber = num2str(i);
    commandExec = strcat(colname,'=','BP_replicate_averaged_by_CTD(:,',colnumber,');');
    evalin('base', commandExec );
end

clear i

% separate data by cruise, then extract data into separate vectors and assign vector
% names 

cruises=unique(Cruise_ID);

for i=1:length(cruises)
    thiscruise=cruises(i);
    cruise_ind=find(BP_replicate_averaged_by_CTD(:,1)==thiscruise);
    subsetname=strcat('BP_replicate_averaged_by_CTD_',num2str(thiscruise));
    for j=1:length(cruise_ind)
        thisrecord=cruise_ind(j);
        commandExec = strcat(subsetname,'(',num2str(j),',',...
            ':)=','BP_replicate_averaged_by_CTD(',num2str(cruise_ind(j)),',:);');
        evalin('base', commandExec );
    end
end

clear i j

% extract data from cruise-specific matrices into separate vectors and then assign vector
% names

for i=1:length(cruises)
    thiscruise=cruises(i);
    cruise_ind=find(BP_replicate_averaged_by_CTD(:,1)==thiscruise);
    subsetname=strcat('BP_replicate_averaged_by_CTD_',num2str(thiscruise));
    for j=1:numcols
        colname = BP_replicate_averaged_by_CTD_headers{1,1}(j);
        colname = char(colname);
        colnumber = num2str(j);
        commandExec = strcat(colname,'_',num2str(thiscruise),'=',subsetname,'(:,',colnumber,');');
        evalin('base', commandExec );
    end
end

clear i j

% now, station data

Station_data_by_station_no = csvread('BLATZ_VICE_BactProd_Station_data_by_station_no.csv',1,0);

numcols = 13;

fid = fopen('BLATZ_VICE_BactProd_Station_data_by_station_no.csv');

Station_data_by_station_no_headers = textscan(fid,'%s',numcols,'delimiter',',');

fclose(fid);

% extract data from matrix into separate vectors and then assign vector
% names

for i=1:numcols
    colname = Station_data_by_station_no_headers{1,1}(i);
    colname = char(colname);
    colnumber = num2str(i);
    commandExec = strcat('Station_data_',colname,'=','Station_data_by_station_no(:,',colnumber,');');
    evalin('base', commandExec );
end

clear i

%% create a depth vs latitude grid and flow in data

depths_2073=unique(Depth_2073);
lats_2073=unique(CTD_station_lat_2073);

Trit_leu_uptake_pmol_L_hr_2073_2D=nan(length(depths_2073),length(lats_2073));

for i=1:length(Trit_leu_uptake_pmol_L_hr_2073)
    if isfinite(Trit_leu_uptake_pmol_L_hr_2073(i))
        thislat=CTD_station_lat_2073(i);
        thisdepth=Depth_2073(i);
        lat_ind=find(lats_2073==thislat);
        depth_ind=find(depths_2073==thisdepth);
        Trit_leu_uptake_pmol_L_hr_2073_2D(depth_ind,lat_ind)=...
                Trit_leu_uptake_pmol_L_hr_2073(i);
    end
end

%% analyze signal:noise by depth

depths=unique(Depth);

signoise_data=nan(100,length(depths));

for i=1:length(depths)
    thisdepth=depths(i);
    depth_ind=find(Depth==thisdepth);
    signoise_data_thisdepth = Signal_to_noise_ratio(depth_ind);
    signoise_data(1:length(signoise_data_thisdepth),i)=signoise_data_thisdepth;
end

% have to get rid of one major outlier

signoise_data(signoise_data>100)=NaN;

nan_ind = find(isnan(signoise_data));
signoise_data_no_blanks = signoise_data;
signoise_data_no_blanks(nan_ind)=[];
signoise_mean = mean(signoise_data_no_blanks(:));

%% some statistical information

subplot(2,2,[1,2]); % uptake rates, all data
boxplot(Trit_leu_uptake_pmol_L_hr,Depth);
xlabel('Depth (m)','FontSize',14);
ylabel(['BP (pmol leu L-1 hr-1)'],'FontSize',14);
xticklabels=[2 8 13 18 24 29 34 40 50 60 70 95 200];
set(gca,'XTick',1:5:depths(length(depths)),'XTickLabel',xticklabels);
title(['Bacterial production rates in the North Atlantic Ocean' char(10)...
    'KN207-1 and KN207-3, all data'],'FontSize',16);

subplot(2,2,3); % live samples and kills, all data
h = errorbar(Depth,Mean_trit_leu_uptake_live_samples_pmol_L_hr,...
    Std_dev_trit_leu_uptake_live_samples_pmol_L_hr,'.','LineWidth',1);
errorbar_tick(h,140);
axis([0 200 0 45]);
xlabel('Depth (m)','FontSize',14);
ylabel(['BP (pmol leu L-1 hr-1)'],'FontSize',14);
hold on
plot(Depth,Trit_leu_uptake_killed_control_pmol_L_hr,'r.');
% set(gca,'XTick',1:5:depths(length(depths)),'XTickLabel',xticklabels);
hleg = legend(['Live samples (mean value; error' char(10) 'bars are standard deviation)'],'Killed control');
set(hleg,'FontSize',9,'box','on');
title(['Live samples vs. killed controls'],'FontSize',14);
hold off;

subplot(2,2,4); % signal to noise
boxplot(signoise_data,depths);
text(58,17,'* mean uptake rate measured in live samples:rate measured in killed control',...
    'HorizontalAlignment','right','FontSize',8);
xlabel('Depth (m)','FontSize',14);
ylabel(['"Signal to noise" ratio*'],'FontSize',14);
xticklabels=[2 8 13 18 24 29 34 40 50 60 70 95 200];
axis([0 62 0 18]);
set(gca,'XTick',1:5:depths(length(depths)),'XTickLabel',xticklabels);
hold on;
plot([0 200],[signoise_mean signoise_mean],'r--');
title(['"Signal to noise" ratio in KN207-1 and KN207-3 bacterial production data'],'FontSize',14);
hold off;

%% convert BP rates to units of C

BP_mmol_C_m3_d_2073 = Trit_leu_uptake_pmol_L_hr_2073*24*1000*(1/10^9).*MW_leu*(1/leu_BPAA)*C_prot*(1/12.01)*ID;
BP_mmol_C_m3_d_2071 = Trit_leu_uptake_pmol_L_hr_2071*24*1000*(1/10^9).*MW_leu*(1/leu_BPAA)*C_prot*(1/12.01)*ID;

BP_mg_C_m3_d_2073 = BP_mmol_C_m3_d_2073*12.01;
BP_mg_C_m3_d_2071 = BP_mmol_C_m3_d_2071*12.01;

BP_mmol_C_m3_d_2073_uncert = Std_dev_trit_leu_uptake_pmol_L_hr_2073*24*1000*(1/10^9).*MW_leu*(1/leu_BPAA)*C_prot*(1/12.01)*ID;
BP_mmol_C_m3_d_2071_uncert = Std_dev_trit_leu_uptake_pmol_L_hr_2071*24*1000*(1/10^9).*MW_leu*(1/leu_BPAA)*C_prot*(1/12.01)*ID;

BP_mg_C_m3_d_2073_uncert = BP_mmol_C_m3_d_2073_uncert*12.01;
BP_mg_C_m3_d_2071_uncert = BP_mmol_C_m3_d_2071_uncert*12.01;

%% KN207-3

figure('Position',[1 scrsz(4)*1/4 scrsz(3)*(10/10) scrsz(4)*4/10])
set(gcf,'PaperPositionMode','auto')

% in mg C/m3/d

% histogram of uptake rates

subplot(2,10,[9:10]);
hist(BP_mg_C_m3_d_2073,30);
xlabel('BP (mg C m-3 d-1)');
ylabel('No. observations');

% contour plot against latitude along cruise track, in mg C m3/d

% first, must interpolate using scatteredInterpolant since data is 
% irregularly spaced; using strictly linear interpolation with contours that connect
% through all data points

ind_exclude=find(isnan(BP_mg_C_m3_d_2073));
BP_mg_C_m3_d_2073_reduced=BP_mg_C_m3_d_2073;
BP_mg_C_m3_d_2073_reduced(ind_exclude)=[];
CTD_station_lat_2073_reduced=CTD_station_lat_2073;
CTD_station_lat_2073_reduced(ind_exclude)=[];
Depth_2073_reduced=Depth_2073;
Depth_2073_reduced(ind_exclude)=[];

F = scatteredInterpolant(CTD_station_lat_2073_reduced,Depth_2073_reduced,...
    BP_mg_C_m3_d_2073_reduced,'natural','linear');

% now, specify some query points and create a grid of data using the
% interpolant

[Xq,Yq] = meshgrid(0:.125:150);

Vq=F(Xq,Yq);

crange=[-0.02 1.8];
numsteps=[-0.02:.025:1.8];
contours=[0 0.5 1.0 1.5 1.8];

subplot(2,10,[1,2,3,4,5,6,7,11,12,13,14,15,16,17])
[c1,h1]=contourf(Xq,-Yq,Vq,numsteps);
hold on
shading flat;
hold on
% [c2,h2]=contour(Xq,-Yq,Vq,contours,'k:');
% clabel(c,h);
% caxis(crange);
xlabel('Degrees north latitude','FontSize',14);
ylabel(['Depth (m)'],'FontSize',14);
axis([CTD_station_lat_2073(1) CTD_station_lat_2073(length(CTD_station_lat_2073))...
    -150 -5]);

h3 = colorbar;
ylabel(h3, 'Bacterial production (mg C m^{-3} d^{-1})')

% overlay CTD data points and color according to whether they've been
% analyzed or not

h4 = plot(CTD_station_lat_2073_reduced,-Depth_2073_reduced,'k+','MarkerSize',5);
h5 = plot(CTD_station_lat_2073(ind_exclude),-Depth_2073(ind_exclude),'ko','MarkerSize',5);
h6 = legend([h4 h5],{'analyzed','no data yet'},'Location','SouthEast','FontSize',12);
set(h6,'FontSize',12,'box','on');
title('Bacterial production along cruise track KN207-3','FontSize',16);
hold off

%% KN207-3

figure('Position',[1 scrsz(4)*1/4 scrsz(3)*(10/10) scrsz(4)*4/10])
set(gcf,'PaperPositionMode','auto')

% in pmol leu/L/hr

% histogram of uptake rates

subplot(2,10,[9:10]);
hist(Trit_leu_uptake_pmol_L_hr_2073,30);
xlabel('BP (pmol leu L-1 h-1)');
ylabel('No. observations');

% contour plot against latitude along cruise track, in pmol leu/L/hr

% first, must interpolate using scatteredInterpolant since data is 
% irregularly spaced; using strictly linear interpolation with contours that connect
% through all data points

ind_exclude=find(isnan(Trit_leu_uptake_pmol_L_hr_2073));
Trit_leu_uptake_pmol_L_hr_2073_reduced=Trit_leu_uptake_pmol_L_hr_2073;
Trit_leu_uptake_pmol_L_hr_2073_reduced(ind_exclude)=[];
CTD_station_lat_2073_reduced=CTD_station_lat_2073;
CTD_station_lat_2073_reduced(ind_exclude)=[];
Depth_2073_reduced=Depth_2073;
Depth_2073_reduced(ind_exclude)=[];

F = scatteredInterpolant(CTD_station_lat_2073_reduced,Depth_2073_reduced,...
    Trit_leu_uptake_pmol_L_hr_2073_reduced,'natural','linear');

% now, specify some query points and create a grid of data using the
% interpolant

[Xq,Yq] = meshgrid(0:.125:150);

Vq=F(Xq,Yq);

crange=[-0.02 1.8]./24/1000/(1/10^9)./MW_leu/(1/leu_BPAA)/C_prot/ID;
numsteps=[-0.02:.025:1.8]./24/1000/(1/10^9)./MW_leu/(1/leu_BPAA)/C_prot/ID;
contours=[0 0.5 1.0 1.5 1.8]./24/1000/(1/10^9)./MW_leu/(1/leu_BPAA)/C_prot/ID;

subplot(2,10,[1,2,3,4,5,6,7,11,12,13,14,15,16,17])
[c1,h1]=contourf(Xq,-Yq,Vq,numsteps);
hold on
shading flat;
hold on
% [c2,h2]=contour(Xq,-Yq,Vq,contours,'k:');
% clabel(c,h);
% caxis(crange);
xlabel('Degrees north latitude','FontSize',14);
ylabel(['Depth (m)'],'FontSize',14);
axis([CTD_station_lat_2073(1) 62 -150 -5]);

h3 = colorbar;
ylabel(h3, 'Bacterial production (pmol leu L^{-1} h^{-1})')

% overlay CTD data points and color according to whether they've been
% analyzed or not

h4 = plot(CTD_station_lat_2073_reduced,-Depth_2073_reduced,'k+','MarkerSize',5);
h5 = plot(CTD_station_lat_2073(ind_exclude),-Depth_2073(ind_exclude),'ko','MarkerSize',5);
h6 = legend([h4 h5],{'analyzed','no data yet'},'Location','SouthEast','FontSize',12);
set(h6,'FontSize',12,'box','on');
title('Bacterial production along cruise track KN207-3','FontSize',16);
hold off

%% KN207-1

figure('Position',[1 scrsz(4)*1/4 scrsz(3)*(10/10) scrsz(4)*4/10])
set(gcf,'PaperPositionMode','auto')

% in mg C/m3/d

% histogram of uptake rates

subplot(2,10,[9:10]);
hist(BP_mg_C_m3_d_2071,30);
xlabel('BP (mg C m-3 d-1)');
ylabel('No. observations');

% contour plot against latitude along cruise track

% first, must interpolate using scatteredInterpolant since data is 
% irregularly spaced; using strictly linear interpolation with contours that connect
% through all data points

ind_exclude=find(isnan(BP_mg_C_m3_d_2071));
BP_mg_C_m3_d_2071_reduced=BP_mg_C_m3_d_2071;
BP_mg_C_m3_d_2071_reduced(ind_exclude)=[];
CTD_station_lat_2071_reduced=CTD_station_lat_2071;
CTD_station_lat_2071_reduced(ind_exclude)=[];
Depth_2071_reduced=Depth_2071;
Depth_2071_reduced(ind_exclude)=[];

F = scatteredInterpolant(CTD_station_lat_2071_reduced,Depth_2071_reduced,...
    BP_mg_C_m3_d_2071_reduced,'natural','linear');

% now, specify some query points and create a grid of data using the
% interpolant

[Xq,Yq] = meshgrid(0:.125:150);

Vq=F(Xq,Yq);

crange=[-0.02 1.8];
numsteps=[-0.02:.025:1.8];
contours=[0 0.5 1.0 1.5 1.8];

subplot(2,10,[1,2,3,4,5,6,7,11,12,13,14,15,16,17])
[c1,h1]=contourf(Xq,-Yq,Vq,numsteps);
hold on
shading flat;
hold on
% [c2,h2]=contour(Xq,-Yq,Vq,contours,'k:');
% clabel(c,h);
% caxis(crange);
xlabel('Degrees north latitude','FontSize',14);
ylabel(['Depth (m)'],'FontSize',14);
axis([CTD_station_lat_2071(length(CTD_station_lat_2071)) CTD_station_lat_2071(1)...
    -150 -5]);

h3 = colorbar;
ylabel(h3, 'Bacterial production (mg C m^{-3} d^{-1})')

% overlay CTD data points and color according to whether they've been
% analyzed or not

h4 = plot(CTD_station_lat_2071_reduced,-Depth_2071_reduced,'k+','MarkerSize',5);
h5 = plot(CTD_station_lat_2071(ind_exclude),-Depth_2071(ind_exclude),'ko','MarkerSize',5);
h6 = legend([h4 h5],{'analyzed','no data yet'},'Location','SouthEast','FontSize',12);
set(h6,'FontSize',12,'box','on');
title('Bacterial production along cruise track KN207-1','FontSize',16);
hold off

%% KN207-1

figure('Position',[1 scrsz(4)*1/4 scrsz(3)*(10/10) scrsz(4)*4/10])
set(gcf,'PaperPositionMode','auto')

% in pmol leu/L/hr

% histogram of uptake rates

subplot(2,10,[9:10]);
hist(Trit_leu_uptake_pmol_L_hr_2071,30);
xlabel('BP (pmol leu L-1 h-1)');
ylabel('No. observations');

% contour plot against latitude along cruise track

% first, must interpolate using scatteredInterpolant since data is 
% irregularly spaced; using strictly linear interpolation with contours that connect
% through all data points

ind_exclude=find(isnan(Trit_leu_uptake_pmol_L_hr_2071));
Trit_leu_uptake_pmol_L_hr_2071_reduced=Trit_leu_uptake_pmol_L_hr_2071;
Trit_leu_uptake_pmol_L_hr_2071_reduced(ind_exclude)=[];
CTD_station_lat_2071_reduced=CTD_station_lat_2071;
CTD_station_lat_2071_reduced(ind_exclude)=[];
Depth_2071_reduced=Depth_2071;
Depth_2071_reduced(ind_exclude)=[];

F = scatteredInterpolant(CTD_station_lat_2071_reduced,Depth_2071_reduced,...
    Trit_leu_uptake_pmol_L_hr_2071_reduced,'natural','linear');

% now, specify some query points and create a grid of data using the
% interpolant

[Xq,Yq] = meshgrid(0:.125:150);

Vq=F(Xq,Yq);

crange=[-0.02 1.8]./24/1000/(1/10^9)./MW_leu/(1/leu_BPAA)/C_prot/ID;
numsteps=[-0.02:.025:1.8]./24/1000/(1/10^9)./MW_leu/(1/leu_BPAA)/C_prot/ID;
contours=[0 0.5 1.0 1.5 1.8]./24/1000/(1/10^9)./MW_leu/(1/leu_BPAA)/C_prot/ID;

subplot(2,10,[1,2,3,4,5,6,7,11,12,13,14,15,16,17])
[c1,h1]=contourf(Xq,-Yq,Vq,numsteps);
hold on
shading flat;
hold on
% [c2,h2]=contour(Xq,-Yq,Vq,contours,'k:');
% clabel(c,h);
% caxis(crange);
xlabel('Degrees north latitude','FontSize',14);
ylabel(['Depth (m)'],'FontSize',14);
axis([CTD_station_lat_2071(length(CTD_station_lat_2071)) CTD_station_lat_2071(1)...
    -150 -5]);

h3 = colorbar;
ylabel(h3, 'Bacterial production (pmol leu L^{-1} h^{-1})')

% overlay CTD data points and color according to whether they've been
% analyzed or not

h4 = plot(CTD_station_lat_2071_reduced,-Depth_2071_reduced,'k+','MarkerSize',5);
h5 = plot(CTD_station_lat_2071(ind_exclude),-Depth_2071(ind_exclude),'ko','MarkerSize',5);
h6 = legend([h4 h5],{'analyzed','no data yet'},'Location','SouthEast','FontSize',12);
set(h6,'FontSize',12,'box','on');
title('Bacterial production along cruise track KN207-1','FontSize',16);
hold off