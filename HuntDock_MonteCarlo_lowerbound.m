% Load Data
close all;
clear all; clc;

load('CO2LAMP_both.mat')
load('SmartRock_both.mat')
load('Sapelovoltage.mat')
load('moonphase_2021.mat')
tides = readtable('C:\Users\lturn\Desktop\all MATLAB files\Sapelo\Sapelo SST\saphdwq2021.csv');% Hunt Dock
%% Combine dates and times to create datetime objects
t = [year(mooninfo2021.e) month(mooninfo2021.Dat) mooninfo2021.VarName1 hour(mooninfo2021.Tim) minute(mooninfo2021.Tim) second(mooninfo2021.Tim)];
d = datetime(t) - hours(5);% subtract 5 hours to convert UTC to ET
%% Find high and low tide times
[pks,locs] = findpeaks(tides.cDepth,tides.DateTimeStamp,'MinPeakDistance',hours(9));
dz = caldays(caldiff(locs))*24; % convert days to hours
hrs = hours(time(caldiff(locs))); % convert difference in time to difference in hours
lowtime = locs(1:end-1) + hours((dz+hrs)./2);% find times when there was a low tide based on half the time between high tides
hightime = locs;% find times when there was a high tide

t = NaT((length(hightime)+length(lowtime)),1);
t(1:2:end) = hightime;
t(2:2:end) = lowtime;
%% Process CO2LAMP data and eliminate unrealistic values
ind = [];
% Replace zeros with data from the next measurement in time
for m = 1:height(co2lampSapelo08172021)
    if co2lampSapelo08172021.CO2ppm(m)==0
        ind = [ind m-1];% record the start of the new measurement cycle
        co2lampSapelo08172021.CO2ppm(m)=co2lampSapelo08172021.CO2ppm(m+1);
    end
end

% Delete values at or below 0
co2lampSapelo08172021.CO2ppm(co2lampSapelo08172021.CO2ppm<=0)=NaN;

Result_tbl = [];
% Create 100 Monte-Carlo style simulations
for j = 1:100
    % Only take last value from each pCO2 measurement cycle
    co2lamp_filtered = co2lampSapelo08172021.CO2ppm(ind(2:end));
    co2lamp_filteredtimes = co2lampSapelo08172021.timestamp(ind(2:end));% now there are 42 data points per day

    for g = 1:length(co2lamp_filtered)
        if ~isnan(co2lamp_filtered(g))
            co2lamp_filtered(g) = round(co2lamp_filtered(g));
            co2lamp_filtered(g) = randi([co2lamp_filtered(g)+20-300 co2lamp_filtered(g)+20]);
        end
    end
    %% Index water temperature data at Lower Duplin River to match timestamps of CO2LAMP data
    DOY = day(co2lamp_filteredtimes,'dayofyear');
    DOY2 = day(tides.DateTimeStamp,'dayofyear');
    q = find(ismember(DOY2,DOY));
    tides2 = tides(q,:);

    % interpolate water temperature to match timestamps of CO2LAMP
    Ttemp = interp1(tides2.DateTimeStamp,tides2.Temp,co2lamp_filteredtimes);
    %% Temp-normalize pCO2
    for idx = 1:length(co2lamp_filtered)
        pCO2_normalized(idx) = co2lamp_filtered(idx)*exp(0.0423*(mean(Ttemp,'omitnan')-Ttemp(idx)));
    end
    pCO2_normalized = pCO2_normalized';
    co2lamp_filtered(co2lamp_filtered>=600)=pCO2_normalized(co2lamp_filtered>=600);%temperature-normalize aquatic pCO2 measurements

    xx = tides2.DateTimeStamp(1):seconds(10):tides2.DateTimeStamp(end);% create more frequent tide times
    yy = spline(datenum(tides2.DateTimeStamp),tides2.cDepth,datenum(xx));% interpolate tide height every ten seconds
    %% Interpolate tide data and eliminate NaN values
    tide = interp1(xx,yy,co2lamp_filteredtimes);% Interpolate tide data to match co2lamp times
    %% Calculating net import/export of CO2
    % Subtract atmospheric CO2 concentration from pCO2
    co2lamp_filtered3 = co2lamp_filtered - 600;
    co2lamp_filtered3(co2lamp_filtered3<=0) = 0;

    count = 1;
    import = nan(size(co2lamp_filteredtimes));
    export = nan(size(co2lamp_filteredtimes));

    for ix = 2:length(tide)
        if tide(ix)>tide(ix-1)
            import(count) = co2lamp_filtered3(ix);
        else
            export(count) = co2lamp_filtered3(ix);
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
    count = 1; % reset counter

    %% Calculate cumulative sums of exports and imports between low tides
    clear import_Gc_cumulative export_Gc_cumulative net_import_Gc export_Gc_cum import_Gc_cum
    count = 1;

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
    net_import_Gc_time = net_import_Gc.*tidecycle';

    Result_tbl(:,j) = [mean(net_import_Gc_time,'omitnan');sum(net_import_Gc_time,'omitnan');max(co2lamp_filtered);min(co2lamp_filtered);nanmean(co2lamp_filtered);std(co2lamp_filtered,'omitnan')/sqrt(length(co2lamp_filtered))];% create a table with 100 variables/columns with all of the results for each run

end

% Calculate the lower-bound of the lower-bound estimates
M = min(Result_tbl,[],2);