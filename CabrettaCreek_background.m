% Load Data
close all;
clear all; clc;

load('co2lamp_detrended2.mat')
load('CO2LAMP_both.mat')
load('SmartRock_both.mat')
load('julyfluxsapelo.mat');
load('moonphase_2021.mat')
tides = readtable('C:\Users\lturn\Desktop\all MATLAB files\Sapelo\Sapelo SST\sapcawq2021.csv');% Cabretta Creek

%% Remove missing value placeholders for flux data
for j= 3:width(flux)
    flux.(j)(flux.(j)<=-9999) = NaN;
end
%% Process CO2LAMP data and eliminate unrealistic values

% Replace zeros with data from the next measurement in time
for m = 1:height(co2lampSapelo08172021)
    if co2lampSapelo08172021.CO2ppm(m)==0
        co2lampSapelo08172021.CO2ppm(m)=co2lampSapelo08172021.CO2ppm(m+1);
    end
end

ind = find(minutes((time(caldiff(co2lampSapelo08172021.timestamp))))>45);
co2lampSapelo08172021 = co2lampSapelo08172021(ind,:);

% Delete values at or below 0
co2lampSapelo08172021.CO2ppm(co2lampSapelo08172021.CO2ppm<=0)=NaN;

DOY = day(co2lamp_filteredtimes,'dayofyear');
DOY2 = day(tides.DateTimeStamp,'dayofyear');
q = find(ismember(DOY2,DOY));
tides2 = tides(q,:);
%% Find high and low tide times
[pks,locs] = findpeaks(tides2.cDepth,tides2.DateTimeStamp,'MinPeakDistance',hours(9));
dz = locs.Day;% convert days to hours
hrs = locs.Hour;% convert difference in time to difference in hours

lowtime = locs + hours((dz+hrs)./2);% find times when there was a low tide based on half the time between high tides
hightime = locs;% find times when there was a high tide

t = NaT((length(hightime)+length(lowtime)),1);
t(1:2:end) = hightime;
t(2:2:end) = lowtime;

Result_tbl = [];

co2lamp_filtered = co2lamp_detrended;
co2lamp_filteredtimes = co2lampSapelo08172021.timestamp;% now there are 22 data points per day
CO2atm = interp1(datenum(flux.TIMESTAMP_END),flux.CO2,datenum(co2lamp_filteredtimes));% interpolate atmospheric CO2 conc.

xx = tides2.DateTimeStamp;
yy = tides2.cDepth;

% Interpolate tide data and eliminate NaN values
tide = interp1(xx,yy,co2lamp_filteredtimes);% Interpolate tide data to match co2lamp times

% Calculating net import/export of CO2
% Subtract atmospheric CO2 concentration from pCO2
co2lamp_filtered3 = co2lamp_filtered - CO2atm;% min

co2lamp_filtered3(co2lamp_filtered3<=0) = 0;

count = 1;
import = nan(size(co2lamp_filteredtimes));
export = nan(size(co2lamp_filteredtimes));

for i = 2:length(tide)
    if tide(i)>tide(i-1)
        import(count) = co2lamp_filtered3(i);
    else
        export(count) = co2lamp_filtered3(i);
    end
    count = count + 1;
end

%% Calculate gC exported per m per second
velocity = abs(diff(tide)./seconds(diff(co2lamp_filteredtimes)));% m/s
export = export(2:end);
import = import(2:end);
% assume air density is 1225 g/m3
export_Gc = (export.*velocity*1225*12)/(28.97*10e6);% gC removed/m2s
import_Gc = (import.*velocity*1225*12)/(28.97*10e6);

%% Calculate cumulative sums of exports and imports between low tides
clear import_Gc_cumulative export_Gc_cumulative net_import_Gc export_Gc_cum import_Gc_cum
count = 1;% reset counter

for i = 1:2:(length(t)-1)
    tf = find(isbetween(co2lamp_filteredtimes,t(i),t(i+1)));
    tf = tf(1:end-1);
    if ~isempty(tf)

        export_Gc_cumulative = cumsum(export_Gc(tf),2,'omitnan');
        export_Gc_cumulative = export_Gc_cumulative(export_Gc_cumulative~=0);

        if ~isempty(export_Gc_cumulative)
            export_Gc_cum(count) = export_Gc_cumulative(end);
        else
            export_Gc_cum(count) = 0;
        end
        import_Gc_cumulative = cumsum(import_Gc(tf),2,'omitnan');
        import_Gc_cumulative = import_Gc_cumulative(import_Gc_cumulative~=0);

        if ~isempty(import_Gc_cumulative)
            import_Gc_cum(count) = import_Gc_cumulative(end);
        else
            import_Gc_cum(count) = 0;
        end

        net_import_Gc(count) = import_Gc_cum(count) - export_Gc_cum(count);
        count = count + 1;
    else
        net_import_Gc(count) = NaN;
        count = count + 1;
    end
end
tidecycle = time(caldiff(t(1:2:end),{'Time'}));
tidecycle = seconds(tidecycle);
net_import_Gc_time = net_import_Gc(1:end-1).*tidecycle';

Result_tbl = [mean(net_import_Gc_time,'omitnan');sum(net_import_Gc_time,'omitnan')];% create a table with results for each run
disp(Result_tbl)

