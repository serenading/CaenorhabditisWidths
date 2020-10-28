close all
clear

%% Script reads width data of different worm strains/species and plots them. 
% DIC images were taken of Day 1 or Day 3 adults on Zeiss Axiophot using 10x objective, 
% and width measurements are taken with FIJI across the worm midbody at three different points.
% Author: @serenading. Oct 2020. 

% set parameters
pixelToMicron = 0.253; % 0.253 um/pixel with 10x objective in large image format (3888 x 2592 pixels). Value given by LMS Microscopy Facility. 
useMax = true; % true or false. If true, script uses the maximum width measurement out of the three for each worm; if false, the average of the three is taken.

%% Get width measurements 
% load  measurement CSV
data = readtable('/Users/sding/OneDrive - Imperial College London/CaenorhabditisWidths/Caenorhabditis_width_measurements.csv');

% get widths dataframe
widths = [data.width_1_pixels,data.width_2_pixels,data.width_3_pixels];

% convert from pixels to microns
widths = widths*pixelToMicron;

% get stats from the three measurements taken for each worm
if useMax
    widths = max(widths,[],2); % use max of the three
else
    widths = mean(widths,2)'; % use average of the three
end

%% Plot double picked Day 1 vs. Day 3 adults
day1LogInd = data.age_days_in_adulthood == 1 & strcmp(data.age_prep_type,'double_pick');
day3LogInd = data.age_days_in_adulthood == 3 & strcmp(data.age_prep_type,'double_pick');
figure; 
subplot(1,2,1)
boxplot(widths(day1LogInd),data.strain(day1LogInd))
ylabel('width (microns)')
title('Day 1 adults')
ylim([0 100])
subplot(1,2,2)
boxplot(widths(day3LogInd),data.strain(day3LogInd))
ylabel('width (microns)')
title('Day 3 adults')
ylim([0 100])
% format
set(findall(gcf,'-property','FontSize'),'FontSize',18)

%% Plot double picked Day 1 rep 1 vs. double picked Day 1 rep 2
rep1LogInd = data.age_days_in_adulthood == 1 & strcmp(data.age_prep_type,'double_pick') & data.date == 20201021;
rep2LogInd = data.age_days_in_adulthood == 1 & strcmp(data.age_prep_type,'double_pick') & data.date == 20201023;
figure; 
subplot(1,2,1)
boxplot(widths(rep1LogInd),data.strain(rep1LogInd))
ylabel('width (microns)')
title('Day 1 adults, rep 1')
ylim([0 100])
subplot(1,2,2)
boxplot(widths(rep2LogInd),data.strain(rep2LogInd))
ylabel('width (microns)')
title('Day 1 adults, rep 2')
ylim([0 100])
% format
set(findall(gcf,'-property','FontSize'),'FontSize',18)
% t-test
strains = unique(cellstr(data.strain));
for strainCtr = 1:numel(strains)
    strain = strains{strainCtr};
    strainLogInd = strcmp(data.strain,strain);
    strainRep1LogInd = strainLogInd & rep1LogInd;
    strainRep2LogInd = strainLogInd & rep2LogInd;
    [h,p] = ttest2(widths(strainRep1LogInd),widths(strainRep2LogInd));
    if h
        disp(['Two-sample t-test returns significant result between the two biological replicates for ' strain '. p value is ' num2str(p) '.'])
    end
end

%% Plot double picked Day 1 vs. bleached Day 1
pickedLogInd = data.age_days_in_adulthood == 1 & strcmp(data.age_prep_type,'double_pick');
bleachedLogInd = data.age_days_in_adulthood == 1 & strcmp(data.age_prep_type,'bleach');
figure; 
subplot(1,2,1)
boxplot(widths(pickedLogInd),data.strain(pickedLogInd))
ylabel('width (microns)')
title('Day 1 adults, double picked')
ylim([0 100])
subplot(1,2,2)
boxplot(widths(bleachedLogInd),data.strain(bleachedLogInd))
ylabel('width (microns)')
title('Day 1 adults, bleached')
ylim([0 100])
% format
set(findall(gcf,'-property','FontSize'),'FontSize',18)
% t-test
N2PickedLogInd = pickedLogInd & strcmp(data.strain,'N2');
[h,p] = ttest2(widths(N2PickedLogInd),widths(bleachedLogInd));
disp(['Two-sample t-test p value is ' num2str(p) ' between double-picked and bleached Day 1 adult N2 worms.'])
