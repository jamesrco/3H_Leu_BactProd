% riBPdata.m: for reading in and working up 3H-leu bacterial production data

% will work up all data as long as it's given to MATLAB in the right
% format, and then export to:
%                  (1) a MATLAB .mat file and
%                  (2) a series of CSV files

% NOTE: non-numerical metadata (notes on each sample or station, etc., get saved
% to a structure within the MATLAB file

% KL 9/12/2011; KL 2/28/2013
% revised JC 11/1/2013

clear all
close all

disp([char(10) 'This script requires the file Quench_curve_fit_params.txt, which you should have generated using plotQuenched.m.']);

% designate a file that will serve as a sink for data output

NameOfFile = 'BLATZ_VICE_BactProd_calcs_20140115.mat';
BP_directory = '/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/BP workups/';
% location of the BP data inputs; also where files wil be dumped at end of
% script

% define and calculate some constants pertaining to the method

Iso_conc_nM = 20; % concentration of 3H-leu in each sample, in nM
Inc_vol_mL = 1; % in mL; 1.0 mL sample used for each incubation
Inc_vol_L = Inc_vol_mL/1000; % in L
Leu_spec_act = 146.5; % specific activity of the 3H used for this work, in Ci/mmol
Ci = 3.7e10; % definition of a Curie, in dps

% read in the LSC data

[num_data txt_data raw_data] = xlsread(strcat(BP_directory,'Van Mooy BP sample master log - KN207-1 and KN207-3 - Ossolinski and Collins 01-14-14.xls'));
data.CPM = cell2mat(raw_data(6:end,9)); % counts per minute, still need to convert to DPM
data.Hnumber = cell2mat(raw_data(6:end,8)); % need this to convert to DPM
data.sInfo = raw_data(6:end,4); % field that contains sample number or notes about the sample if it's a standard or blank

% also need the sample number since that will be used to merge with the
% information in the sample inventory; some samples dont have a number
% (e.g., blanks and quenching standards) so have to assign as NaN

data.sNumber(1:size(data.sInfo,1),1) = NaN;
for i = 1:size(data.sInfo,1);
    if isnumeric(cell2mat(data.sInfo(i)));
       data.sNumber(i,1) = cell2mat(data.sInfo(i));
    else
        data.sNumber(i,1) = NaN;
    end
end
clear i

% now, read in data from the BP sample inventory; will need to merge this
% with our LSC data since the sample inventory has the details for each BP vial
% JC transcribed this data from the two BP field notebooks in 10/2013
% will have a mix of CTD data and data from some on-deck mesocosm
% experiments

% the inventory spreadsheet has two tabs, first is a sample inventory
% containing information about each sample; second is a station inventory
% containing information about each station

% first, read in data from the sample inventory

[num_sampinv txt_sampinv raw_sampinv] = xlsread(strcat(BP_directory,'Sample inventory for BP workups - KN207-1,KN207-3.xlsx'),'Sample ID inventory');

% clean things up a bit

raw_sampinv=raw_sampinv(5:end,:);
ind_badsamples = find(cellfun(@(x) any(isnan(x)), raw_sampinv(:,2)));
raw_sampinv(ind_badsamples,:)=[];
num_sampinv(ind_badsamples,:)=[];

disp([char(10) 'The following sample IDs will be excluded from data processing because they contain no data. Recommend checking data input to see whether this should be the case.' char(10)]);
disp(ind_badsamples);

% now, load data

sampinv.sNumber = cell2mat(raw_sampinv(:,1)); % sample number from the inventory
                                          % this is the field on which we'll 
                                          % eventually join the two arrays
sampinv.CruiseName = raw_sampinv(:,2);
sampinv.CruiseID(1:size(sampinv.sNumber,1),1) = NaN;
ind = find(cellfun(@(x) strcmp(x,'KN207-1'), sampinv.CruiseName));
sampinv.CruiseID(ind)=2071;
ind = find(cellfun(@(x) strcmp(x,'KN207-3'), sampinv.CruiseName));
sampinv.CruiseID(ind)=2073;
sampinv.Sam_type = raw_sampinv(:,3); % sample type, CTD cast or incubation experiment
sampinv.CTDStation = cell2mat(raw_sampinv(:,4)); % a number (for experiments, will be NaN)
                                             % will be used to join sample
                                             % inventory on station
                                             % inventory
sampinv.Depth = cell2mat(raw_sampinv(:,5)); % in meters (for experiments, will be NaN)
sampinv.Exper_name = raw_sampinv(:,6); % experiment name; for CTD stations will be NaN
sampinv.Exper_timept = cell2mat(raw_sampinv(:,7)); % sampling timepoint for experiments 
sampinv.Exper_treat = raw_sampinv(:,8); % experimental treatment applied to this sample
sampinv.Exper_treat_carboy = cell2mat(raw_sampinv(:,9)); % experiment carboy number
sampinv.LiveKilled = raw_sampinv(:,10); % live = nothing added, dead = killed with 5% TCA
sampinv.LiveKilledindex(1:size(sampinv.sNumber,1),1) = NaN; % live=1 killed=2
ind = find(cellfun(@(x) strcmp(x,'l'), sampinv.LiveKilled));
sampinv.LiveKilledindex(ind)=1;
ind = find(cellfun(@(x) strcmp(x,'k'), sampinv.LiveKilled));
sampinv.LiveKilledindex(ind)=2;
sampinv.Notes = raw_sampinv(:,11);

% convert experiment names to numbers to make data processing easier
sampinv.Exper_no(1:size(sampinv.sNumber,1),1) = NaN;
ind_exp1=find(cellfun(@(x) strcmp(x,'Process_Station_1_mesocosms'), sampinv.Exper_name));
sampinv.Exper_no(ind_exp1) = 1;
ind_exp2=find(cellfun(@(x) strcmp(x,'Process_Station_2_mesocosms'), sampinv.Exper_name));
sampinv.Exper_no(ind_exp2) = 2;
ind_exp3=find(cellfun(@(x) strcmp(x,'H2O2_addition_mesocosms'), sampinv.Exper_name));
sampinv.Exper_no(ind_exp3) = 3;
ind_exp4=find(cellfun(@(x) strcmp(x,'Cocco_1_mesocosms'), sampinv.Exper_name));
sampinv.Exper_no(ind_exp4) = 4;
ind_exp5=find(cellfun(@(x) strcmp(x,'Cocco_1_bullseye_H2O2_addition_mesocosms'), sampinv.Exper_name));
sampinv.Exper_no(ind_exp5) = 5;
ind_exp6=find(cellfun(@(x) strcmp(x,'Ehv_mesocosms'), sampinv.Exper_name));
sampinv.Exper_no(ind_exp6) = 6;
ind_exp7=find(cellfun(@(x) strcmp(x,'Hi_lo_light_mesocosms'), sampinv.Exper_name));
sampinv.Exper_no(ind_exp7) = 7;

% convert experiment treatments to numbers to make data processing easier
sampinv.Exper_treat_no(1:size(sampinv.sNumber,1),1) = NaN;
ind_tmt1=find(cellfun(@(x) strcmp(x,'control'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt1) = 1;
ind_tmt2=find(cellfun(@(x) strcmp(x,'plus_nutrients'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt2) = 2;
ind_tmt3=find(cellfun(@(x) strcmp(x,'plus_myriocin'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt3) = 3;
ind_tmt4=find(cellfun(@(x) strcmp(x,'temp'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt4) = 4;
ind_tmt5=find(cellfun(@(x) strcmp(x,'temp_plus_nutrients'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt5) = 5;
ind_tmt6=find(cellfun(@(x) strcmp(x,'plus_5_uM_H2O2'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt6) = 6;
ind_tmt7=find(cellfun(@(x) strcmp(x,'plus_30_uM_H2O2'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt7) = 7;
ind_tmt8=find(cellfun(@(x) strcmp(x,'plus_10_uM_H2O2'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt8) = 8;
ind_tmt9=find(cellfun(@(x) strcmp(x,'plus_50_uM_H2O2'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt9) = 9;
ind_tmt10=find(cellfun(@(x) strcmp(x,'plus_EhV1'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt10) = 10;
ind_tmt11=find(cellfun(@(x) strcmp(x,'plus_EhV1_plus_myriocin'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt11) = 11;
ind_tmt12=find(cellfun(@(x) strcmp(x,'low_light_control'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt12) = 12;
ind_tmt13=find(cellfun(@(x) strcmp(x,'low_light_plus_EhV1'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt13) = 13;
ind_tmt14=find(cellfun(@(x) strcmp(x,'high_light_control'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt14) = 14;
ind_tmt15=find(cellfun(@(x) strcmp(x,'high_light_plus_EhV1'), sampinv.Exper_treat));
sampinv.Exper_treat_no(ind_tmt15) = 15;

% now, read in data from the station inventory

[num_stainv txt_stainv  raw_stainv ] = xlsread(strcat(BP_directory,'Sample inventory for BP workups - KN207-1,KN207-3.xlsx'),'Incu durations, temps, notes');

stainv.CruiseName = raw_stainv(4:end,1);
stainv.CruiseID(1:size(stainv.CruiseName,1),1) = NaN;
ind = find(cellfun(@(x) strcmp(x,'KN207-1'), stainv.CruiseName));
stainv.CruiseID(ind)=2071;
ind = find(cellfun(@(x) strcmp(x,'KN207-3'), stainv.CruiseName));
stainv.CruiseID(ind)=2073;
stainv.Sam_type = raw_stainv(4:end,2); % sample type, CTD cast or incubation experiment
stainv.CTDStation = cell2mat(raw_stainv(4:end,3)); % a number (for experiments, will be NaN)
                                             % will be used to join sample
                                             % inventory on station
                                             % inventory
stainv.Exper_name = raw_stainv(4:end,4); % experiment name; for CTD stations will be NaN
stainv.Exper_timept = cell2mat(raw_stainv(4:end,5)); % sampling timepoint for experiments 
stainv.IncTime_start = cell2mat(raw_stainv(4:end,6)); % start of the incubation time
stainv.IncTime_end = cell2mat(raw_stainv(4:end,7)); % end of the incubation time
stainv.CTDincutemp = cell2mat(raw_stainv(4:end,8)); % in deg C, incubation temp
                                                    % for CTD station samples
stainv.Exper_baseincutemp = cell2mat(raw_stainv(4:end,9)); % in deg C, incubation temp
                                                    % for experimental samples - all except temp treatment
stainv.Exper_temptreatementincutemp = cell2mat(raw_stainv(4:end,10));
% in deg C, incubation temp for experimental samples to which temp
% treatment was applied
stainv.Station_Notes = raw_stainv(4:end,11);

% convert experiment names to numbers to make data processing easier
stainv.Exper_no(1:size(stainv.CruiseID,1),1) = NaN;
ind_exp1_stainv=find(cellfun(@(x) strcmp(x,'Process_Station_1_mesocosms'), stainv.Exper_name));
stainv.Exper_no(ind_exp1_stainv) = 1;
ind_exp2_stainv=find(cellfun(@(x) strcmp(x,'Process_Station_2_mesocosms'), stainv.Exper_name));
stainv.Exper_no(ind_exp2_stainv) = 2;
ind_exp3_stainv=find(cellfun(@(x) strcmp(x,'H2O2_addition_mesocosms'), stainv.Exper_name));
stainv.Exper_no(ind_exp3_stainv) = 3;
ind_exp4_stainv=find(cellfun(@(x) strcmp(x,'Cocco_1_mesocosms'), stainv.Exper_name));
stainv.Exper_no(ind_exp4_stainv) = 4;
ind_exp5_stainv=find(cellfun(@(x) strcmp(x,'Cocco_1_bullseye_H2O2_addition_mesocosms'), stainv.Exper_name));
stainv.Exper_no(ind_exp5_stainv) = 5;
ind_exp6_stainv=find(cellfun(@(x) strcmp(x,'Ehv_mesocosms'), stainv.Exper_name));
stainv.Exper_no(ind_exp6_stainv) = 6;
ind_exp7_stainv=find(cellfun(@(x) strcmp(x,'Hi_lo_light_mesocosms'), stainv.Exper_name));
stainv.Exper_no(ind_exp7_stainv) = 7;

% read in CTD station data

[num_KN2071data txt_KN2071data  raw_KN2071data ] = xlsread(strcat(BP_directory,'KN207-1_shipcasts.xlsx'));
[num_KN2073data txt_KN2073data  raw_KN2073data ] = xlsread(strcat(BP_directory,'KN207-3_shipcasts.xlsx'));

Ship_data_2071=cell2mat(raw_KN2071data(9:end,[1 5 6 7]));
Ship_data_2073=cell2mat(raw_KN2073data(11:end,[1 5 6 7]));

Ship_data = nan(length(Ship_data_2071)+length(Ship_data_2073),5);

Ship_data(1:length(Ship_data_2071),1)=2071;
Ship_data(length(Ship_data_2071)+1:end,1)=2073;

Ship_data(1:length(Ship_data_2071),2:5)=Ship_data_2071;
Ship_data(length(Ship_data_2071)+1:end,2:5)=Ship_data_2073;

% now, calculate the duration of each incubation

stainv.IncTime_dur(1:size(stainv.IncTime_start,1),1) = NaN;

IncTime_start=stainv.IncTime_start;
IncTime_end=stainv.IncTime_end;

for i = 1:size(IncTime_start,1)
    t1 = IncTime_start(i,1);
    t2 = IncTime_end(i,1);

    stainv.IncTime_dur(i,1) = (t2-t1)*24*60; % duration in minutes
    
    clear t1 t2
end

clear i IncTime_start IncTime_end

% now, to get values in DPM rather than CPM,
% run the following after getting the numbers from plotQuenched.m

p = load('Quench_curve_fit_params.txt');
x = data.Hnumber;
data.efficiency = p(1).*(x.^3) - p(2).*(x.^2) - p(3).*x + p(4) ;
data.DPM = data.CPM./data.efficiency; clear x

% next, let's join the LSC data matrix and sample inventory on each other

[sNum ia ib]=intersect(data.sNumber,sampinv.sNumber);

% now, put the data in a results matrix; non-numerical data can go later into a "notes" structure

BP_results = nan(length(sampinv.sNumber),24);

BP_results(:,1)=sampinv.sNumber;
BP_results(:,2)=sampinv.CruiseID;
BP_results(:,3)=sampinv.CTDStation;
BP_results(:,4)=sampinv.Depth;
BP_results(:,5)=sampinv.Exper_no;
BP_results(:,6)=sampinv.Exper_timept;
BP_results(:,7)=sampinv.Exper_treat_no;
BP_results(:,8)=sampinv.Exper_treat_carboy;
BP_results(:,9)=sampinv.LiveKilledindex;

% feed in LSC data

BP_results(ib,10)=data.DPM(ia);
BP_results(ib,11)=data.CPM(ia);
BP_results(ib,12)=data.Hnumber(ia);
BP_results(ib,13)=data.efficiency(ia);

% put numerical station inventory information into one matrix

Station_data = nan(length(stainv.CruiseID),13);

Station_data(:,1)=stainv.CruiseID;
Station_data(:,2)=stainv.CTDStation;
Station_data(:,3)=stainv.Exper_no;
Station_data(:,4)=stainv.Exper_timept;
Station_data(:,5)=stainv.IncTime_start;
Station_data(:,6)=stainv.IncTime_end;
Station_data(:,7)=stainv.IncTime_dur;
Station_data(:,8)=stainv.CTDincutemp;
Station_data(:,9)=stainv.Exper_baseincutemp;
Station_data(:,10)=stainv.Exper_temptreatementincutemp;

% get some information about the dataset

cruises = unique(BP_results(:,2));
CTD_stations = unique(BP_results(:,3));
CTD_stations(isnan(CTD_stations))=[];
depths = unique(BP_results(:,4));
depths(isnan(depths))=[];
experiments = unique(BP_results(:,5));
experiments(isnan(experiments))=[];
exper_timepoints= unique(BP_results(:,6));
exper_timepoints(isnan(exper_timepoints))=[];
exper_treatments= unique(BP_results(:,7));
exper_treatments(isnan(exper_treatments))=[];

% use a for loop to match information from Ship_data (CTD cast time and
% lat/long) with CTD station number in the Station_data matrix, and then
% store that information in some new columns in the Station_data

for i=1:length(cruises)
    thiscruise=cruises(i);
        for j=1:length(CTD_stations)
            thisstation=CTD_stations(j);
            Station_data_ind=find(Station_data(:,1)==thiscruise & Station_data(:,2)==thisstation);
            Ship_data_ind=find(Ship_data(:,1)==thiscruise & Ship_data(:,2)==thisstation);
            if isfinite(Station_data_ind)
            Station_data(Station_data_ind,11)=Ship_data(Ship_data_ind,3);
            Station_data(Station_data_ind,12)=Ship_data(Ship_data_ind,4);
            Station_data(Station_data_ind,13)=Ship_data(Ship_data_ind,5);
            end
        end
end

clear i j

% now, use a series of for loops to feed in relevant station data for each sample ID
% into the right spot in the "BP_results" data matrix

for i=1:length(cruises)
    
    thiscruise=cruises(i);
    
    ind_thiscruise_BP_results=find(BP_results(:,2)==thiscruise);
    BP_results_thiscruise=BP_results(ind_thiscruise_BP_results,:);
    ind_thiscruise_Station_data=find(Station_data(:,1)==thiscruise);
    Station_data_thiscruise=Station_data(ind_thiscruise_Station_data,:);
    CTDstasthiscruise=unique(Station_data(ind_thiscruise_Station_data,2));
    CTDstasthiscruise(isnan(CTDstasthiscruise))=[];
    
    for j=1:length(CTDstasthiscruise)
        thisstation=CTDstasthiscruise(j);
        ind_thisstation_BP_results=find(BP_results_thiscruise(:,3)==thisstation);
        ind_thisstation_Station_data=find(Station_data_thiscruise(:,2)==thisstation);
        BP_results_thiscruise(ind_thisstation_BP_results,14)=...
            Station_data_thiscruise(ind_thisstation_Station_data,5);
        BP_results_thiscruise(ind_thisstation_BP_results,15)=...
            Station_data_thiscruise(ind_thisstation_Station_data,6);
        BP_results_thiscruise(ind_thisstation_BP_results,16)=...
            Station_data_thiscruise(ind_thisstation_Station_data,7);
        BP_results_thiscruise(ind_thisstation_BP_results,17)=...
            Station_data_thiscruise(ind_thisstation_Station_data,8);
        BP_results_thiscruise(ind_thisstation_BP_results,18)=...
            Station_data_thiscruise(ind_thisstation_Station_data,9);
        BP_results_thiscruise(ind_thisstation_BP_results,19)=...
            Station_data_thiscruise(ind_thisstation_Station_data,10);
         BP_results_thiscruise(ind_thisstation_BP_results,22)=...
            Station_data_thiscruise(ind_thisstation_Station_data,11);
        BP_results_thiscruise(ind_thisstation_BP_results,23)=...
            Station_data_thiscruise(ind_thisstation_Station_data,12);
        BP_results_thiscruise(ind_thisstation_BP_results,24)=...
            Station_data_thiscruise(ind_thisstation_Station_data,13);       
        
        [sNum ia ib]=intersect(BP_results_thiscruise(:,1),BP_results(:,1));
        
        BP_results(ib,14:19)=BP_results_thiscruise(ia,14:19);
        BP_results(ib,22:24)=BP_results_thiscruise(ia,22:24);

    end
    
    Experimentsthiscruise=unique(Station_data(ind_thiscruise_Station_data,3));
    Experimentsthiscruise(isnan(Experimentsthiscruise))=[];
    
    for k=1:length(Experimentsthiscruise)
        thisexperiment=Experimentsthiscruise(k);

        ind_thisexperiment_BP_results=find(BP_results_thiscruise(:,5)==thisexperiment);
        ind_thisexperiment_Station_data=find(Station_data_thiscruise(:,3)==thisexperiment);
        BP_results_thisexperiment=BP_results_thiscruise(ind_thisexperiment_BP_results,:);
        Station_data_thisexperiment=Station_data_thiscruise(ind_thisexperiment_Station_data,:);

        Timeptsthisexpt=unique(Station_data_thiscruise(ind_thisexperiment_Station_data,4));
        Timeptsthisexpt(isnan(Timeptsthisexpt))=[];  
        
        for l=1:length(Timeptsthisexpt)
            
            thistimept=Timeptsthisexpt(l);
            ind_thistimept_BP_results=find(BP_results_thisexperiment(:,6)==thistimept);
            ind_thistimept_Station_data=find(Station_data_thisexperiment(:,4)==thistimept);
 
            BP_results_thisexperiment(ind_thistimept_BP_results,14)=...
            Station_data_thisexperiment(ind_thistimept_Station_data,5);
            BP_results_thisexperiment(ind_thistimept_BP_results,15)=...
            Station_data_thisexperiment(ind_thistimept_Station_data,6);
            BP_results_thisexperiment(ind_thistimept_BP_results,16)=...
            Station_data_thisexperiment(ind_thistimept_Station_data,7);
            BP_results_thisexperiment(ind_thistimept_BP_results,17)=...
            Station_data_thisexperiment(ind_thistimept_Station_data,8);
            BP_results_thisexperiment(ind_thistimept_BP_results,18)=...
            Station_data_thisexperiment(ind_thistimept_Station_data,9);
            BP_results_thisexperiment(ind_thistimept_BP_results,19)=...
            Station_data_thisexperiment(ind_thistimept_Station_data,10);
        
            [sNum ia ib]=intersect(BP_results_thisexperiment(:,1),BP_results_thiscruise(:,1));
        
            BP_results_thiscruise(ib,14:19)=BP_results_thisexperiment(ia,14:19);
        end
        
        [sNum ia ib]=intersect(BP_results_thiscruise(:,1),BP_results(:,1));
        
        BP_results(ib,14:19)=BP_results_thiscruise(ia,14:19);
        BP_results(ib,22:24)=BP_results_thiscruise(ia,22:24);

    end
        
end

clear i j k l

% calculate uptake rates

Trit_leu_uptake_dpm_L_min=BP_results(:,10)*(1/Inc_vol_L).*(1./BP_results(:,16)); % in dpm/L/min

Trit_leu_uptake_dpm_L_hr=Trit_leu_uptake_dpm_L_min*60;

Trit_leu_uptake_pmol_L_hr=Trit_leu_uptake_dpm_L_hr*(1/Ci)*(1/60)*(1/Leu_spec_act)*10^9; % in pmol/L/hr

; % in pmol/L/hr

BP_results(:,20)=Trit_leu_uptake_dpm_L_hr;
BP_results(:,21)=Trit_leu_uptake_pmol_L_hr;

% average the 3 live reps, calculate 3H-leu uptake rates, and attach
% station-specific data to each CTD station depth and experimental
% treatment

% first for CTD stations

% Will yield a results matrix "BP_results_replicate_averaged_by_CTD" 
% with the following structure:
%
% Cruise_ID CTD_station_no Depth Mean_DPM_live_samples Standard_dev_live_samples
% DPM_of_killed_control Incu_duration Mean_trit_leu_uptake_live_samples_dpm_L_hr
% Std_dev_trit_leu_uptake_live_samples_dpm_l_hr Trit_leu_uptake_killed_control_dpm_L_hr
% Trit_leu_uptake_dpm_L_hr Std_dev_trit_leu_uptake_dpm_l_hr
% Mean_trit_leu_uptake_live_samples_pmol_L_hr Std_dev_trit_leu_uptake_live_samples_pmol_L_hr
% Trit_leu_uptake_killed_control_pmol_L_hr Trit_leu_uptake_pmol_L_hr Std_dev_trit_leu_uptake_pmol_L_hr
% Signal_to_noise_ratio CTD_sample_incutemp

BP_results_replicate_averaged_by_CTD = nan(1,22);

i=1;

for j=1:length(cruises)
    thiscruise=cruises(j);
    for k=1:length(CTD_stations)
        thisstation=CTD_stations(k);
        for l=1:length(depths)
            thisdepth=depths(l);
            killedcontrol_ind=find(BP_results(:,2)==thiscruise & BP_results(:,3)==thisstation...
                & BP_results(:,4)==thisdepth & BP_results(:,9)==2);
            livereps_ind=find(BP_results(:,2)==thiscruise & BP_results(:,3)==thisstation...
                & BP_results(:,4)==thisdepth & BP_results(:,9)==1);
            if isfinite(livereps_ind)
                mean_DPM_live=mean(BP_results(livereps_ind,10));
                stddev_DPM_live=std(BP_results(livereps_ind,10));
                DPM_killed=BP_results(killedcontrol_ind,10);
                Mean_trit_leu_uptake_live_samples_dpm_L_hr=mean(BP_results(livereps_ind,20));
                Std_dev_trit_leu_uptake_live_samples_dpm_l_hr=std(BP_results(livereps_ind,20));
                Trit_leu_uptake_killed_control_dpm_L_hr=BP_results(killedcontrol_ind,20);
                Trit_leu_uptake_dpm_L_hr=Mean_trit_leu_uptake_live_samples_dpm_L_hr-Trit_leu_uptake_killed_control_dpm_L_hr;
                Std_dev_trit_leu_uptake_dpm_l_hr=Std_dev_trit_leu_uptake_live_samples_dpm_l_hr;
                Mean_trit_leu_uptake_live_samples_pmol_L_hr=mean(BP_results(livereps_ind,21));
                Std_dev_trit_leu_uptake_live_samples_pmol_L_hr=std(BP_results(livereps_ind,21));
                Trit_leu_uptake_killed_control_pmol_L_hr=BP_results(killedcontrol_ind,21);
                Trit_leu_uptake_pmol_L_hr=Mean_trit_leu_uptake_live_samples_pmol_L_hr-Trit_leu_uptake_killed_control_pmol_L_hr;
                Std_dev_trit_leu_uptake_pmol_L_hr=Std_dev_trit_leu_uptake_live_samples_pmol_L_hr;
                Signal_to_noise_ratio=Mean_trit_leu_uptake_live_samples_pmol_L_hr/Trit_leu_uptake_killed_control_pmol_L_hr;
                BP_results_replicate_averaged_by_CTD(i,1:22)=[thiscruise thisstation ...
                    thisdepth mean_DPM_live stddev_DPM_live DPM_killed ...
                    BP_results(killedcontrol_ind,16) ...
                    Mean_trit_leu_uptake_live_samples_dpm_L_hr ...
                    Std_dev_trit_leu_uptake_live_samples_dpm_l_hr ...
                    Trit_leu_uptake_killed_control_dpm_L_hr ...
                    Trit_leu_uptake_dpm_L_hr ...
                    Std_dev_trit_leu_uptake_dpm_l_hr ...
                    Mean_trit_leu_uptake_live_samples_pmol_L_hr ...
                    Std_dev_trit_leu_uptake_live_samples_pmol_L_hr ...
                    Trit_leu_uptake_killed_control_pmol_L_hr ...
                    Trit_leu_uptake_pmol_L_hr ...
                    Std_dev_trit_leu_uptake_pmol_L_hr ...
                    Signal_to_noise_ratio ...
                    BP_results(killedcontrol_ind,17) ...
                    BP_results(killedcontrol_ind,22) ...
                    BP_results(killedcontrol_ind,23) ...
                    BP_results(killedcontrol_ind,24)];
                i=i+1;
            else
            end
        end
    end
end

clear i j k l

% then for experiments

BP_results_replicate_averaged_by_experiment = nan(1,24);

% Will yield a results matrix "BP_results_replicate_averaged_by_experiment" 
% with the following structure:
%
% Cruise_ID Experiment_no Experiment_time_point Exper_treatment_no Mean_DPM_live_samples
% Standard_dev_live_samples Mean_DPM_killed_controls
% Standard_dev_killed_controls Incu_duration Mean_trit_leu_uptake_live_samples_dpm_L_hr
% Std_dev_trit_leu_uptake_live_samples_dpm_l_hr Mean_trit_leu_uptake_killed_controls_dpm_L_hr
% Std_dev_trit_leu_uptake_killed_controls_dpm_L_hr 
% Trit_leu_uptake_dpm_L_hr Std_dev_trit_leu_uptake_dpm_l_hr Mean_trit_leu_uptake_live_samples_pmol_L_hr
% Std_dev_trit_leu_uptake_live_samples_pmol_L_hr Mean_trit_leu_uptake_killed_controls_pmol_L_hr
% Std_dev_trit_leu_uptake_killed_controls_pmol_L_hr 
% Trit_leu_uptake_pmol_L_hr Std_dev_trit_leu_uptake_pmol_L_hr
% Signal_to_noise_ratio Experimental_sample_incutemp Experimental_temp_treatment_incu_temp

i=1;

for j=1:length(cruises)
    thiscruise=cruises(j);
    for k=1:length(experiments)
        thisexpt=experiments(k);
        for l=1:length(exper_timepoints)
            thistimept=exper_timepoints(l);
            for m=1:length(exper_treatments)
                thistreatment=exper_treatments(m);
                killedcontrols_ind=find(BP_results(:,2)==thiscruise & BP_results(:,5)==thisexpt...
                & BP_results(:,6)==thistimept & BP_results(:,7)==thistreatment & BP_results(:,9)==2);
                livereps_ind=find(BP_results(:,2)==thiscruise & BP_results(:,5)==thisexpt...
                & BP_results(:,6)==thistimept & BP_results(:,7)==thistreatment & BP_results(:,9)==1);
            if isfinite(livereps_ind)
                mean_DPM_live=mean(BP_results(livereps_ind,10));
                stddev_DPM_live=std(BP_results(livereps_ind,10));
                mean_DPM_killed=mean(BP_results(killedcontrols_ind,10));
                stddev_DPM_killed=std(BP_results(killedcontrols_ind,10));
                Mean_trit_leu_uptake_live_samples_dpm_L_hr=mean(BP_results(livereps_ind,20));
                Std_dev_trit_leu_uptake_live_samples_dpm_l_hr=std(BP_results(livereps_ind,20));
                Mean_trit_leu_uptake_killed_controls_dpm_L_hr=mean(BP_results(killedcontrols_ind,20));
                Std_dev_trit_leu_uptake_killed_controls_dpm_L_hr=std(BP_results(killedcontrols_ind,20));
                Trit_leu_uptake_dpm_L_hr=Mean_trit_leu_uptake_live_samples_dpm_L_hr-Mean_trit_leu_uptake_killed_controls_dpm_L_hr;
                Std_dev_trit_leu_uptake_dpm_l_hr=sqrt(Std_dev_trit_leu_uptake_live_samples_dpm_l_hr^2 + ...
                    Std_dev_trit_leu_uptake_killed_controls_dpm_L_hr^2);       
                Mean_trit_leu_uptake_live_samples_pmol_L_hr=mean(BP_results(livereps_ind,21));
                Std_dev_trit_leu_uptake_live_samples_pmol_L_hr=std(BP_results(livereps_ind,21));
                Mean_trit_leu_uptake_killed_controls_pmol_L_hr=mean(BP_results(killedcontrols_ind,21));
                Std_dev_trit_leu_uptake_killed_controls_pmol_L_hr=std(BP_results(killedcontrols_ind,21));
                Trit_leu_uptake_pmol_L_hr=Mean_trit_leu_uptake_live_samples_pmol_L_hr-Mean_trit_leu_uptake_killed_controls_pmol_L_hr;
                Std_dev_trit_leu_uptake_pmol_L_hr=sqrt(Std_dev_trit_leu_uptake_live_samples_pmol_L_hr^2 + ...
                    Std_dev_trit_leu_uptake_killed_controls_pmol_L_hr^2);
                Signal_to_noise_ratio=Mean_trit_leu_uptake_live_samples_pmol_L_hr/Mean_trit_leu_uptake_killed_controls_pmol_L_hr;
                BP_results_replicate_averaged_by_experiment(i,1:24)=[thiscruise thisexpt ...
                    thistimept thistreatment mean_DPM_live stddev_DPM_live mean_DPM_killed ...
                    stddev_DPM_killed BP_results(killedcontrols_ind(1),16) ...
                    Mean_trit_leu_uptake_live_samples_dpm_L_hr ...
                    Std_dev_trit_leu_uptake_live_samples_dpm_l_hr ...
                    Mean_trit_leu_uptake_killed_controls_dpm_L_hr ...
                    Std_dev_trit_leu_uptake_killed_controls_dpm_L_hr ...
                    Trit_leu_uptake_dpm_L_hr ...
                    Std_dev_trit_leu_uptake_dpm_l_hr ...
                    Mean_trit_leu_uptake_live_samples_pmol_L_hr ...
                    Std_dev_trit_leu_uptake_live_samples_pmol_L_hr ...
                    Mean_trit_leu_uptake_killed_controls_pmol_L_hr ...
                    Std_dev_trit_leu_uptake_killed_controls_pmol_L_hr ...
                    Trit_leu_uptake_pmol_L_hr ...
                    Std_dev_trit_leu_uptake_pmol_L_hr ...
                    Signal_to_noise_ratio ...
                    BP_results(killedcontrols_ind(1),18) ...
                    BP_results(killedcontrols_ind(1),19)];
                i=i+1;
            else
            end
            end
        end
    end
end
    
clear i j k l m

% now, clean things up
    
clear txt_stainv BP_results_thiscruise BP_results_thisexperiment CTDstasthiscruise Experimentsthiscruise ...
    Station_data_thiscruise Station_data_thisexperiment Timeptsthisexpt cruises data i ia ib ind ...
    ind_thiscruise_BP_results ind_thiscruise_Station_data ind_thisexperiment_BP_results ...
    ind_thisexperiment_Station_data ind_thisstation_BP_results ind_thisstation_Station_data ...
    ind_thistimept_BP_results ind_thistimept_Station_data ind_tmt1 ind_tmt10 ind_tmt11 ind_tmt1 ...
    ind_tmt13 ind_tmt14 ind_tmt15 ind_tmt2 ind_tmt3 ind_tmt4 ind_tmt5 ind_tmt6 ind_tmt7 ind_tmt8 ind_tmt9 ...
    j k l num_data num_sampinv num_stainv p raw_data raw_sampinv raw_stainv sNum thiscruise ...
    thisexperiment thisstation thistimept txt_data txt_sampinv ...
    ind_exp1 ind_exp1_stainv ind_exp2 ind_exp2_stainv ind_exp3 ind_exp3_stainv ind_exp4 ...
    ind_exp4_stainv ind_exp5 ind_exp5_stainv ind_exp6 ind_exp6_stainv ind_exp7 ind_exp7_stainv ind_tmt12 ...
    DPM_killed killedcontrol_ind killedcontrols_ind livereps_ind mean_DPM_killed mean_DPM_live ...
    stddev_DPM_killed stddev_DPM_live thisdepth thisexpt thistreatment Trit_leu_uptake_dpm_L_hr ...
    Trit_leu_uptake_dpm_L_min Trit_leu_uptake_pmol_L_hr ind_badsamples  txt_KN2073data ...
    txt_KN2071data raw_KN2073data raw_KN2071data outid num_KN2073data num_KN2071data ...
    headers_Station_data_by_station_no headers_BP_results_replicate_averaged_by_experiment ...
    headers_BP_results_replicate_averaged_by_CTD headers_BP_results_by_sample_ID export_precision ...
    ans Trit_leu_uptake_killed_control_pmol_L_hr Trit_leu_uptake_killed_control_dpm_L_hr ...
    Std_dev_trit_leu_uptake_pmol_L_hr Std_dev_trit_leu_uptake_live_samples_pmol_L_hr ...
    Std_dev_trit_leu_uptake_live_samples_dpm_l_hr Std_dev_trit_leu_uptake_killed_controls_pmol_L_hr ...
    Std_dev_trit_leu_uptake_killed_controls_dpm_L_hr Std_dev_trit_leu_uptake_dpm_l_hr Station_data_ind ...
    Signal_to_noise_ratio Ship_data_ind Ship_data_2073 Ship_data_2071 Ship_data ...
    Mean_trit_leu_uptake_live_samples_pmol_L_hr Mean_trit_leu_uptake_live_samples_dpm_L_hr ...
    Mean_trit_leu_uptake_killed_controls_pmol_L_hr Mean_trit_leu_uptake_killed_controls_dpm_L_hr ...
    experiments depths exper_timepoints exper_treatments Leu_spec_act Iso_conc_nM Inc_vol_mL ...
    Inc_vol_L Ci CTD_stations

% rename our two "big" matrices to something more descriptive

BP_results_by_sample_ID=BP_results;
Station_data_by_station_no=Station_data;

% create a new cell array with the non-numerical data

BP_results_by_sample_ID_with_metadata.Sample_number=BP_results(:,1);
BP_results_by_sample_ID_with_metadata.Cruise_name=sampinv.CruiseName;
BP_results_by_sample_ID_with_metadata.Cruise_ID=BP_results(:,2);
BP_results_by_sample_ID_with_metadata.Sample_type=sampinv.Sam_type;
BP_results_by_sample_ID_with_metadata.CTD_station_no=BP_results(:,3);
BP_results_by_sample_ID_with_metadata.Depth=BP_results(:,4);
BP_results_by_sample_ID_with_metadata.Experiment_name=sampinv.Exper_name;
BP_results_by_sample_ID_with_metadata.Experiment_no=BP_results(:,5);
BP_results_by_sample_ID_with_metadata.Experiment_time_point=BP_results(:,6);
BP_results_by_sample_ID_with_metadata.Experiment_treatment_descr=sampinv.Exper_treat;
BP_results_by_sample_ID_with_metadata.Exper_treatment_no=BP_results(:,7);
BP_results_by_sample_ID_with_metadata.Exper_treatment_carboy_no=BP_results(:,8);
BP_results_by_sample_ID_with_metadata.Live_or_killed=sampinv.LiveKilled;
BP_results_by_sample_ID_with_metadata.Live_or_killed_index=BP_results(:,9);
BP_results_by_sample_ID_with_metadata.DPM=BP_results(:,10);
BP_results_by_sample_ID_with_metadata.CPM=BP_results(:,11);
BP_results_by_sample_ID_with_metadata.H_number=BP_results(:,12);
BP_results_by_sample_ID_with_metadata.LSC_efficiency=BP_results(:,13);
BP_results_by_sample_ID_with_metadata.Incu_start_time=BP_results(:,14);
BP_results_by_sample_ID_with_metadata.Incu_end_time=BP_results(:,15);
BP_results_by_sample_ID_with_metadata.Incu_duration=BP_results(:,16);
BP_results_by_sample_ID_with_metadata.Trit_leu_uptake_dpm_L_hr=BP_results(:,20);
BP_results_by_sample_ID_with_metadata.Trit_leu_uptake_pmol_L_hr=BP_results(:,21);
BP_results_by_sample_ID_with_metadata.CTD_sample_incutemp=BP_results(:,17);
BP_results_by_sample_ID_with_metadata.Experimental_sample_incutemp=BP_results(:,18);
BP_results_by_sample_ID_with_metadata.Experimental_temp_treatment_incu_temp=BP_results(:,19);
BP_results_by_sample_ID_with_metadata.CTD_station_time=BP_results(:,22);
BP_results_by_sample_ID_with_metadata.CTD_station_lat=BP_results(:,23);
BP_results_by_sample_ID_with_metadata.CTD_station_long=BP_results(:,24);
BP_results_by_sample_ID_with_metadata.Sample_notes=sampinv.Notes;

% do the same for station information

Station_info_with_metadata.Cruise_name=stainv.CruiseName;
Station_info_with_metadata.Cruise_ID=Station_data(:,1);
Station_info_with_metadata.Sample_type=stainv.Sam_type;
Station_info_with_metadata.CTD_station_no=Station_data(:,2);
Station_info_with_metadata.Experiment_name=stainv.Exper_name;
Station_info_with_metadata.Experiment_no=Station_data(:,3);
Station_info_with_metadata.Incu_start_time=Station_data(:,5);
Station_info_with_metadata.Incu_end_time=Station_data(:,6);
Station_info_with_metadata.Incu_duration=Station_data(:,7);
Station_info_with_metadata.CTD_sample_incutemp=Station_data(:,8);
Station_info_with_metadata.Experimental_sample_incutemp=Station_data(:,9);
Station_info_with_metadata.Experimental_temp_treatment_incu_temp=Station_data(:,10);
Station_info_with_metadata.CTD_station_time=Station_data(:,11);
Station_info_with_metadata.CTD_station_lat=Station_data(:,12);
Station_info_with_metadata.CTD_station_long=Station_data(:,13);
Station_info_with_metadata.Station_notes=stainv.Station_Notes;

% clean up some more stuff

clear BP_results Station_data sampinv stainv

% save to a MATLAB file

save(strcat(BP_directory,NameOfFile));

% export to a series of CSV files

export_precision=10; % necessary otherwise timestamps won't be written correctly

headers_BP_results_by_sample_ID = ['Sample_number,Cruise_ID,CTD_station_no,Depth_m,Experiment_no,Experiment_time_point,Exper_treat_no,Exper_treatment_carboy_no,Live_or_killed_index,DPM,CPM,H_number,LSC_efficiency,Incu_start_time,Incu_end_time,Incu_duration,CTD_sample_incutemp,Experimental_sample_incutemp,Experimental_temp_treatment_incu_temp,Trit_leu_uptake_dpm_L_hr,Trit_leu_uptake_pmol_L_hr,CTD_station_time,CTD_station_lat,CTD_station_long'];
outid = fopen(strcat(BP_directory,'BLATZ_VICE_BactProd_results_by_sample_ID.csv'), 'w+');
fprintf(outid, '%s', headers_BP_results_by_sample_ID);
fclose(outid);
dlmwrite (strcat(BP_directory,'BLATZ_VICE_BactProd_results_by_sample_ID.csv'),BP_results_by_sample_ID,'roffset',1,'-append','precision',export_precision);

headers_Station_data_by_station_no = ['Cruise_ID,CTD_station_no,Experiment_no,Experiment_time_point,Incu_start_time,Incu_end_time,Incu_duration,CTD_sample_incutemp,Experimental_sample_incutemp,Experimental_temp_treatment_incu_temp,CTD_station_time,CTD_station_lat,CTD_station_long'];
outid = fopen(strcat(BP_directory,'BLATZ_VICE_BactProd_Station_data_by_station_no.csv'), 'w+');
fprintf(outid, '%s', headers_Station_data_by_station_no);
fclose(outid);
dlmwrite (strcat(BP_directory,'BLATZ_VICE_BactProd_Station_data_by_station_no.csv'),Station_data_by_station_no,'roffset',1,'-append','precision',export_precision);

headers_BP_results_replicate_averaged_by_CTD = ['Cruise_ID,CTD_station_no,Depth,Mean_DPM_live_samples,Standard_dev_live_samples,DPM_of_killed_control,Incu_duration,Mean_trit_leu_uptake_live_samples_dpm_L_hr,Std_dev_trit_leu_uptake_live_samples_dpm_l_hr,Trit_leu_uptake_killed_control_dpm_L_hr,Trit_leu_uptake_dpm_L_hr,Std_dev_trit_leu_uptake_dpm_l_hr,Mean_trit_leu_uptake_live_samples_pmol_L_hr,Std_dev_trit_leu_uptake_live_samples_pmol_L_hr,Trit_leu_uptake_killed_control_pmol_L_hr,Trit_leu_uptake_pmol_L_hr,Std_dev_trit_leu_uptake_pmol_L_hr,Signal_to_noise_ratio,CTD_sample_incutemp,CTD_station_time,CTD_station_lat,CTD_station_long'];
outid = fopen(strcat(BP_directory,'BLATZ_VICE_BactProd_results_replicate_averaged_by_CTD.csv'), 'w+');
fprintf(outid, '%s', headers_BP_results_replicate_averaged_by_CTD);
fclose(outid);
dlmwrite (strcat(BP_directory,'BLATZ_VICE_BactProd_results_replicate_averaged_by_CTD.csv'),BP_results_replicate_averaged_by_CTD,'roffset',1,'-append','precision',export_precision);

headers_BP_results_replicate_averaged_by_experiment = ['Cruise_ID,Experiment_no,Experiment_time_point,Exper_treatment_no,Mean_DPM_live_samples,Standard_dev_live_samples,Mean_DPM_killed_controls,Standard_dev_killed_controls,Incu_duration,Mean_trit_leu_uptake_live_samples_dpm_L_hr,Std_dev_trit_leu_uptake_live_samples_dpm_l_hr,Mean_trit_leu_uptake_killed_controls_dpm_L_hr,Std_dev_trit_leu_uptake_killed_controls_dpm_L_hr,Trit_leu_uptake_dpm_L_hr,Std_dev_trit_leu_uptake_dpm_l_hr,Mean_trit_leu_uptake_live_samples_pmol_L_hr,Std_dev_trit_leu_uptake_live_samples_pmol_L_hr,Mean_trit_leu_uptake_killed_controls_pmol_L_hr,Std_dev_trit_leu_uptake_killed_controls_pmol_L_hr,Trit_leu_uptake_pmol_L_hr,Std_dev_trit_leu_uptake_pmol_L_hr,Signal_to_noise_ratio,Experimental_sample_incutemp,Experimental_temp_treatment_incu_temp'];
outid = fopen(strcat(BP_directory,'BLATZ_VICE_BactProd_results_replicate_averaged_by_experiment.csv'), 'w+');
fprintf(outid, '%s', headers_BP_results_replicate_averaged_by_experiment);
fclose(outid);
dlmwrite (strcat(BP_directory,'BLATZ_VICE_BactProd_results_replicate_averaged_by_experiment.csv'),BP_results_replicate_averaged_by_experiment,'roffset',1,'-append','precision',export_precision);

disp([char(10) 'The following files have been generated:' char(10)]);
disp('BLATZ_VICE_BactProd_results_by_sample_ID.csv');
disp('BLATZ_VICE_BactProd_Station_data_by_station_no.csv');
disp('BLATZ_VICE_BactProd_results_replicate_averaged_by_CTD.csv');
disp('BLATZ_VICE_BactProd_results_replicate_averaged_by_experiment.csv');
disp(NameOfFile);
disp([char(10) 'These files can be found in:']);
disp(BP_directory);

% clean up some more stuff

clear outid headers_Station_data_by_station_no headers_BP_results_replicate_averaged_by_experiment ...
    headers_BP_results_replicate_averaged_by_CTD headers_BP_results_by_sample_ID export_precision ...
    BP_directory NameOfFile

