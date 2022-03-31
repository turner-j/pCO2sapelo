% Load Data
close all;
clear all; clc;

load('CO2LAMP_both.mat')
load('SmartRock_both.mat')
load('Sapelovoltage.mat')
load('moonphase_2021.mat')
tides = readtable('C:\Users\lturn\Desktop\all MATLAB files\Sapelo\Sapelo SST\sapcawq2021.csv');% Cabretta Creek
latlon = readtable('C:\Users\lturn\Desktop\all MATLAB files\Sapelo\Sapelo SST\sampling_stations.csv');
%% Combine dates and times to create datetime objects
t = [year(mooninfo2021.e) month(mooninfo2021.Dat) mooninfo2021.VarName1 hour(mooninfo2021.Tim) minute(mooninfo2021.Tim) second(mooninfo2021.Tim)];
d = datetime(t) - hours(5);% subtract 5 hours to convert UTC to ET
%% Process CO2LAMP data and eliminate unrealistic values
ind = [];
% Replace zeros with data from the next measurement in time
for i=1:height(co2lampSapelo08172021)
    if co2lampSapelo08172021.CO2ppm(i)==0
        ind = [ind i-1];% record the start of the new measurement cycle
        co2lampSapelo08172021.CO2ppm(i)=co2lampSapelo08172021.CO2ppm(i+1);
    end
end

% Delete values at or below 0
co2lampSapelo08172021.CO2ppm(co2lampSapelo08172021.CO2ppm<=0)=NaN;
% co2lampSapelo08172021.CO2ppm(co2lampSapelo08172021.CO2ppm>10000)=NaN;

% Only take last value from each pCO2 measurement cycle
co2lamp_filtered = co2lampSapelo08172021.CO2ppm(ind(2:end));
co2lamp_filteredtimes = co2lampSapelo08172021.timestamp(ind(2:end));% now there are 42 data points per day

%% Index water temperature data at Lower Duplin River to match timestamps of CO2LAMP data
DOY = day(co2lamp_filteredtimes,'dayofyear');
DOY2 = day(tides.DateTimeStamp,'dayofyear');
ind = find(ismember(DOY2,DOY));
tides = tides(ind,:);

Ttemp = interp1(tides.DateTimeStamp,tides.Temp,co2lamp_filteredtimes);% interpolate water temperature to match timestamps of CO2LAMP
Tsal = interp1(tides.DateTimeStamp,tides.Sal,co2lamp_filteredtimes);% repeat with salinity
TpH = interp1(tides.DateTimeStamp,tides.pH,co2lamp_filteredtimes);% pH
Tturb = interp1(tides.DateTimeStamp,tides.Turb,co2lamp_filteredtimes);% turbidity (NTU)
TDO = interp1(tides.DateTimeStamp,tides.DO_mgl,co2lamp_filteredtimes);% DO (mg/L)
TSpCond = interp1(tides.DateTimeStamp,tides.SpCond,co2lamp_filteredtimes);% SpCond (mS/cm)
   
%% Temp-normalize pCO2
for i = 1:length(co2lamp_filtered)
    pCO2_normalized(i) = co2lamp_filtered(i)*exp(0.0423*(mean(Ttemp,'omitnan')-Ttemp(i)));
end
pCO2_normalized = pCO2_normalized';
co2lamp_filtered(co2lamp_filtered>=600)=pCO2_normalized(co2lamp_filtered>=600);%temperature-normalize aquatic pCO2 measurements

max(co2lamp_filtered)
min(co2lamp_filtered)
nanmean(co2lamp_filtered)
stderror = std(co2lamp_filtered,'omitnan')/sqrt(length(co2lamp_filtered))

xx = tides.DateTimeStamp(1):seconds(10):tides.DateTimeStamp(end);% create more frequent tide times
yy = spline(datenum(tides.DateTimeStamp),tides.cDepth,datenum(xx));% interpolate tide height every ten seconds
%% Interpolate tide data and eliminate NaN values
tide = interp1(xx,yy,co2lamp_filteredtimes);% Interpolate tide data to match co2lamp times

% bin average according to tide
co2lamp_filtered2 = co2lamp_filtered;
co2lamp_filtered2(co2lamp_filtered2<600) = NaN;% delete values not taken underwater

tide = tide(~isnan(co2lamp_filtered2),:);
co2lamp_filtered2(isnan(co2lamp_filtered2))=[];

% Calculate some basic statistics
max(co2lamp_filtered2)
min(co2lamp_filtered2)
nanmean(co2lamp_filtered2)
stderror = std(co2lamp_filtered2)/sqrt(length(co2lamp_filtered2))

%% Create boxplot/boxchart of pCO2 at different tide heights
binEdges = min(tide):.2:max(tide);
grouptide = discretize(tide,binEdges);
b = boxchart(grouptide,co2lamp_filtered2)
% Calculate mean of each bin
% meanpCO2 = groupsummary(co2lamp_filtered2,grouptide,'mean');
bns = (binEdges(1:end-1) + binEdges(2:end))./2;
y1 = -828.3.*(bns.^2)+2602.*(bns)+819.4;
hold on
plot(y1,'-k')
ylabel('{\itp}CO_2')
xlabel('Tide Height (m)')
formatSpec = '%.2f';
text(0.2,5500,[{'R^2 ='} num2str(0.86,formatSpec)],'FontSize',12)
caption = sprintf('y = %.0f*x^2 +%.0f*x + %.0f',-828.3,2602,819.4);
text(0.2,6500, caption, 'FontSize', 12, 'Color', 'k');
% xlim([0 14])
% xticks(1:13)
% xticklabels([0 round(bns(2:end),1)])
ylim([500 7000])
box on
set(gca,'FontSize',17)

b.BoxFaceColor = [0.5333    0.6039    0.8902];
b.BoxFaceAlpha = 0.1;
b.LineWidth = 1.5;
b.WhiskerLineColor =  [0.2039    0.2941    0.6588];
b.MarkerColor = [0.5333    0.6039    0.8902];

ind = find((tide>1.29)&(tide<1.49));
range(co2lamp_filtered2(ind))
std(co2lamp_filtered2(ind))

ind = find((tide>=2.09));
range(co2lamp_filtered2(ind))
std(co2lamp_filtered2(ind))

% figure()
% plot(binEdges,meanpCO2,'.')
%% Interpolate co2lamp to match interpolated tide data
pCO2 = interp1(co2lamp_filteredtimes,co2lamp_filtered,xx);
tide = interp1(xx,yy,co2lamp_filteredtimes);

pCO2(isnan(pCO2))=0;
yy(isnan(yy))=0;

%% Analyze grainger causality of links between pCO2 and NEE of CO2 due to tide
vq4 = interp1(xx,yy,co2lamp_filteredtimes); % tide height
phase = interp1(d,mooninfo2021.Phase,co2lamp_filteredtimes); % moon phase
% switch up the order to see if the results change
tbl = table(phase./10,vq4,Tsal,TSpCond,Tturb,TpH,TDO,Ttemp,co2lamp_filtered./100); % create a table with the response as the last variable
tbl.Properties.VariableNames = {'Moon Phase','Tide Height','Salinity','Specific Conductivity','Turbidity','pH','Dissolved Oxygen','Temp','pCO2'};
%% Fit VAR models, with lags ranging from 1 to 4, to the series. Initialize each fit by specifying the first four observations. Store the Akaike information criteria (AIC) of the fits.
T = size(tbl,1); % Total sample size
numseries = size(tbl,2);
numlags = (1:(numseries+1))';
nummdls = numel(numlags);

% Partition time base.
maxp = max(numlags); % Maximum number of required presample responses
idxpre = 1:maxp;
idxest = (maxp + 1):T;

% Preallocation
EstMdl(nummdls) = varm(numseries,0);
aic = zeros(nummdls,1);

% Fit VAR models to data.
Y0 = tbl{idxpre,:}; % Presample
Y = tbl{idxest,:};  % Estimation sample
for j = 1:numel(numlags)
    Mdl = varm(numseries,numlags(j));
    Mdl.SeriesNames = tbl.Properties.VariableNames;
    EstMdl(j) = estimate(Mdl,Y,'Y0',Y0);
    results = summarize(EstMdl(j));
    aic(j) = results.AIC;
end

[~,bestidx] = min(aic);
% Select the model that yields the best fit
BestMdl = EstMdl(bestidx);
%% Plot observations and fitted values on same graph
resid = infer(BestMdl,Y);

yhat = tbl.pCO2(10:end) - resid(:,3);

co2lampcropped = co2lamp_filtered(10:end);
[R,P] = corrcoef(100*yhat(~isnan(yhat)),(co2lampcropped(~isnan(yhat))))

%%
[h,Summary] = gctest(BestMdl)

Mdl = fitrensemble(tbl,'pCO2');
imp = predictorImportance(Mdl);
figure
bar(imp)
title('Predictor Importance Estimates')
ylabel('Estimates')
xlabel('Predictors')
ax = gca;
ax.XTickLabel = Mdl.PredictorNames;
%% Subplots
fig = figure();
colorsnew = colormap(parula(10));

subplot(4,1,1)
hold on
plot(co2lamp_filteredtimes,co2lamp_filtered,'color',colorsnew(1,:),'LineWidth',2)
% plot(co2lamp_filteredtimes(10:end),100*yhat,'r:','LineWidth',2)
% plot(co2lampSapelo08172021.CO2ppm)
hold off
ylabel('{\it p}CO_2')% (ppm)
set(gca,'FontSize',17)
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
set(gca,'xticklabel',[])
xticks(hightime)

subplot(4,1,2)
plot(data8172021.timestamp,data8172021.ECcorrected,'color',colorsnew(3,:),'LineWidth',2)
ylabel('salinity')% (\muS/cm), is there a limit on the salinity sensor?
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
set(gca,'FontSize',17)
set(gca,'xticklabel',[])
xticks(hightime)

voltage = (data8172021.turbidity2.*3.3)/4096;

subplot(4,1,3)
plot(data8172021.timestamp,voltage,'color',colorsnew(4,:),'LineWidth',2)
ylabel('turbidity')%  (JTU)
set(gca,'FontSize',17)
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
xticks(hightime)
set(gca,'xticklabel',[])

subplot(4,1,4)
plot(tides.DateTimeStamp,tides.cDepth,'.','Color',colorsnew(5,:))
hold on
plot(xx,yy,'Color',colorsnew(6,:))
ylabel('height')%(m)
set(gca,'FontSize',17)
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
xtickformat('D')
xticks(hightime)
xlabel('Day of Year')

%% Calculating net import/export of CO2
% Subtract atmospheric CO2 concentration from pCO2
co2lamp_filtered = co2lamp_filtered - 600;
co2lamp_filtered(co2lamp_filtered<=0) = 0;

count = 1;
import = nan(size(co2lamp_filteredtimes));
export = nan(size(co2lamp_filteredtimes));

for i = 2:length(tide)
    if tide(i)>tide(i-1)
        import(count) = co2lamp_filtered(i);
    else
        export(count) = co2lamp_filtered(i);
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
clear import_Gc_cumulative export_Gc_cumulative net_export_Gc export_Gc_cum import_Gc_cum
count = 1;

t = NaT((length(hightime)+length(lowtime)),1);
t(1:2:end) = hightime;
t(2:2:end) = lowtime;

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
mean(net_import_Gc_time,'omitnan')
sum(net_import_Gc_time,'omitnan')
%% Create map to display tide gauge locations
lat = latlon.Latitude;
lon = latlon.Longitude;
gb = geobubble(lat,-lon,[],categorical(latlon.StationName),'Basemap','landcover');
% gb.BubbleWidthRange = 25;
gb.MapLayout = 'maximized';
% gb.ZoomLevel = 14;
geobasemap colorterrain
% addToolbarMapButton(gb,"basemap")