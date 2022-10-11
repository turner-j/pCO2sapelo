% Load Data
close all;
clear ; clc;

load('CO2LAMP_both.mat')
load('SmartRock_both.mat')
load('Sapelovoltage.mat')
load('moonphase_2021.mat')
load('julyfluxsapelo.mat');
% tide data from https://tidesandcurrents.noaa.gov/noaatidepredictions.html?id=8675622&units=standard&bdate=20210817&edate=20210917&timezone=LST/LDT&clock=12hour&datum=MLLW&interval=hilo&action=dailychart
% Import water temp:
tides = readtable('E:\SanDiskSecureAccess\all MATLAB files\Sapelo\Sapelo SST\saphdwq2021.csv');% Hunt Dock
%% Combine dates and times to create datetime objects
t = [year(mooninfo2021.e) month(mooninfo2021.Dat) mooninfo2021.VarName1 hour(mooninfo2021.Tim) minute(mooninfo2021.Tim) second(mooninfo2021.Tim)];
d = datetime(t) - hours(5);% subtract 5 hours to convert UTC to ET

% Remove missing value placeholders for flux data
for j= 3:width(flux)
    flux.(j)(flux.(j)<=-9999) = NaN;
end
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
% Replace zeros with data from the next measurement in time
for i=1:height(co2lampSapelo08172021)
    if co2lampSapelo08172021.CO2ppm(i)==0
        co2lampSapelo08172021.CO2ppm(i)=co2lampSapelo08172021.CO2ppm(i+1);
    end
end

ind = find(minutes((time(caldiff(co2lampSapelo08172021.timestamp))))>45);
co2lampSapelo08172021 = co2lampSapelo08172021(ind,:);

% Delete values at or below 0
co2lampSapelo08172021.CO2ppm(co2lampSapelo08172021.CO2ppm<=0)=NaN;

co2lamp_filtered = co2lampSapelo08172021.CO2ppm;
co2lamp_filteredtimes = co2lampSapelo08172021.timestamp;
CO2atm = interp1(datenum(flux.TIMESTAMP_END),flux.CO2,datenum(co2lamp_filteredtimes));% interpolate atmospheric CO2 conc.

% % Calculate how many points per day
[~,~,ix] = unique(day(co2lamp_filteredtimes,'dayofyear'));
C = median(accumarray(ix,1).');

%% Detrend pCO2
figure()
plot(movmin(co2lamp_filtered(~isnan(co2lamp_filtered)),C),'.')
h = lsline;
test = polyfit(get(h,'xdata'),get(h,'ydata'),1);
slopem = test(1);intercpt = test(2);

figure()
plot(co2lamp_filteredtimes(~isnan(co2lamp_filtered)),movmin(co2lamp_filtered(~isnan(co2lamp_filtered)),C),'o','LineWidth',2,'Color',[0.6549    0.8078    0.9098]) 
hold on

for i = 1:length(co2lamp_filteredtimes)
    co2trend(i) = (slopem.*i)+intercpt;
end

plot(co2lamp_filteredtimes,co2trend,'--','LineWidth',2,'Color',[0.6549    0.8078    0.9098])
caption = sprintf('y = %.1f*n + %.0f',slopem,intercpt);
text(co2lamp_filteredtimes(round(.5*length(co2lamp_filteredtimes))),560, caption, 'FontSize', 15, 'Color',[0.6549    0.8078    0.9098]);
ylabel('{\itp}CO_2 moving min.')
set(gca,'FontSize',15)

[R,P]=corrcoef(co2trend(~isnan(co2lamp_filtered)),co2lamp_filtered(~isnan(co2lamp_filtered)))
co2lamp_detrended = co2lamp_filtered - co2trend'+420; % detrend all data
co2lamp_filtered = co2lamp_detrended;

plot(co2lamp_filteredtimes(~isnan(co2lamp_filtered)),movmin(co2lamp_filtered(~isnan(co2lamp_filtered)),C),'^','Color',[0.0588    0.3098    0.4784])
plot(co2lamp_filteredtimes,ones(size(co2lamp_filteredtimes)).*420,'--','LineWidth',2,'Color',[0.0588    0.3098    0.4784])
caption = sprintf('y = %.0f',420);
text(co2lamp_filteredtimes(round(.66*length(co2lamp_filteredtimes))),450, caption, 'FontSize', 15, 'Color',[0.0588    0.3098    0.4784]);
hold off
legend('original','best fit','detrended','best fit','Location','northwest')
%% Index water temperature data at Lower Duplin River to match timestamps of CO2LAMP data
DOY = day(co2lamp_filteredtimes,'dayofyear');
DOY2 = day(tides.DateTimeStamp,'dayofyear');
ind = find(ismember(DOY2,DOY));
tides = tides(ind,:);

% interpolate tide data to match timestamps of CO2LAMP
Tsal = interp1(tides.DateTimeStamp,tides.Sal,co2lamp_filteredtimes);% repeat with salinity
Tturb = interp1(tides.DateTimeStamp,tides.Turb,co2lamp_filteredtimes);% turbidity (NTU)
TSpCond = interp1(tides.DateTimeStamp,tides.SpCond,co2lamp_filteredtimes);% SpCond (mS/cm)
Ttemp = interp1(tides.DateTimeStamp,tides.Temp,co2lamp_filteredtimes);% water temp (deg C)
Atemp = interp1(flux.TIMESTAMP_END,flux.TA_EP,co2lamp_filteredtimes);% water temp (deg C)

[R,P]=corrcoef(Ttemp(~isnan(co2lamp_filtered)),co2lamp_filtered(~isnan(co2lamp_filtered)))

% min(tides.Temp)
% max(tides.Temp)

% max(co2lamp_filtered)
% min(co2lamp_filtered)
% nanmean(co2lamp_filtered)
% stderror = std(co2lamp_filtered,'omitnan')/sqrt(length(co2lamp_filtered))

xx = tides.DateTimeStamp;
yy = tides.cDepth;

%% Interpolate tide data and eliminate NaN values
tide = interp1(xx,yy,co2lamp_filteredtimes);% Interpolate tide data to match co2lamp times

% remove background CO2 using flux tower data
co2lamp_filtered2 = co2lamp_filtered-CO2atm;

tide2 = tide(co2lamp_filtered2>0);
co2lamp_filtered2(co2lamp_filtered2<=0)=[];% delete values not taken underwater
co2lamp_filtered2(isnan(co2lamp_filtered2))=[];

% Calculate some basic statistics
max(co2lamp_filtered2)
min(co2lamp_filtered2)
nanmean(co2lamp_filtered2)
stderror = std(~isnan(co2lamp_filtered2))/sqrt(length(~isnan(co2lamp_filtered2)))

%% Create boxplot/boxchart of pCO2 at different tide heights
binEdges = min(tide2):.2:max(tide2);
grouptide = discretize(tide2,binEdges);

% Calculate mean of each bin
bns = (binEdges(1:end-1) + binEdges(2:end))./2;

%Calculate mean of each bin
medianpCO2 = groupsummary(co2lamp_filtered2,grouptide,'median');
mdl = fitlm(binEdges,medianpCO2,'Intercept',false);% Calculate best fit line
y1 = mdl.Coefficients.Estimate(1).*bns;

figure()
hold on
b = boxchart(grouptide,co2lamp_filtered2);
plot(y1,'-k')
ylabel('{\itp}CO_2')
xlabel('Tide Height (m)')
caption = sprintf('R^2 = %.2f',mdl.Rsquared.Ordinary);
text(11,4200,caption,'FontSize',12, 'Color', 'k')
caption = sprintf('y = %.0f*x',mdl.Coefficients.Estimate(1));
text(9,5000, caption, 'FontSize', 12, 'Color', 'k');
xlim([0 16])
xticks(1:15)
xticklabels([round(bns,1)])
ylim([0 6000])
box on
set(gca,'FontSize',17)

b.BoxFaceColor = [0.5333    0.6039    0.8902];
b.BoxFaceAlpha = 0.1;
b.LineWidth = 1.5;
b.WhiskerLineColor =  [0.2039    0.2941    0.6588];
b.MarkerColor = [0.5333    0.6039    0.8902];
hold off

for i = 1:15
    dtst = co2lamp_filtered2((ismember(grouptide,i))==1);
    co2lampgroup(i) = sum(isoutlier(dtst,'quartiles')); % count number of outliers for each box
end

%% Interpolate co2lamp to match interpolated tide data
pCO2 = interp1(co2lamp_filteredtimes,co2lamp_filtered,xx);
tide = interp1(xx,yy,co2lamp_filteredtimes);

pCO2(isnan(pCO2))=0;
yy(isnan(yy))=0;
%% Multiple linear regression
vq4 = interp1(xx,yy,co2lamp_filteredtimes); % tide height
phase = interp1(d,mooninfo2021.Phase,co2lamp_filteredtimes); % moon phase

tbl = table(vq4,Tsal./10,Tturb./10,Ttemp./10,co2lamp_filtered./100); % create a table with the response as the last variable
tbl.Properties.VariableNames = {'TideHeight','Salinity','Turbidity','Temp','pCO2'};
tbl2 = tbl;
% Withhold half of the data to test model fit
% tbl(1:2:end,:) = [];% odd-numbered observations only
%% Fit VAR models, with lags, to the series. Initialize each fit by specifying the first four observations. Store the Akaike information criteria (AIC) of the fits.
% Odd-numbered observations only
% T = size(tbl,1); % Total sample size
% numseries = size(tbl,2);
% numlags = (1:(numseries+1))';
% nummdls = numel(numlags); 
% 
% % Partition time base.
% maxp = max(numlags); % Maximum number of required presample responses
% idxpre = 1:maxp;
% idxest = (maxp + 1):T;
% 
% % Preallocation
% EstMdl(nummdls) = varm(numseries,0);
% aic = zeros(nummdls,1);
% 
% % Fit VAR models to data.
% Y0 = tbl{idxpre,:}; % Presample
% Y = tbl{idxest,:};  % Estimation sample
% for j = 1:numel(numlags)
%     Mdl = varm(numseries,numlags(j));
%     Mdl.SeriesNames = tbl.Properties.VariableNames;
%     EstMdl(j) = estimate(Mdl,Y,'Y0',Y0);
%     results = summarize(EstMdl(j));
%     aic(j) = results.AIC;
% end
% 
% [~,bestidx] = min(aic);
% % Select the model that yields the best fit
% BestMdl = EstMdl(bestidx);

%% Calculate model fit to other half of data
tbl3 = tbl2;
% tbl2(2:2:end,:) = [];% even-numbered observations only
% resid = infer(BestMdl,table2array(tbl2));% calculate the residuals when using withheld data
% 
% % Moonphase_yhat = tbl2.MoonPhase((length(tbl2.pCO2)-length(resid)+1):end)
% % - resid(:,1);% estimate moon phase
% Tideht_yhat = tbl2.TideHeight((length(tbl2.pCO2)-length(resid)+1):end) - resid(:,1);% estimate tide ht
% Salinity_yhat = tbl2.Salinity((length(tbl2.pCO2)-length(resid)+1):end) - resid(:,2);% estimate salinity
% Turbidity_yhat = tbl2.Turbidity((length(tbl2.pCO2)-length(resid)+1):end) - resid(:,3);% estimate turb
% Temp_yhat = tbl2.Temp((length(tbl2.pCO2)-length(resid)+1):end) - resid(:,4);% estimate temp
% yhat = tbl2.pCO2((length(tbl2.pCO2)-length(resid)+1):end) - resid(:,5);% estimate pCO2
% 
% co2lampcropped = 100.*tbl2.pCO2((length(tbl2.pCO2)-length(resid)+1):end);
% [R,P] = corrcoef(100*yhat(~isnan(yhat)),(co2lampcropped(~isnan(yhat)))) % How well VAR fits pCO2 obs

%% Now with full dataset
% Fit VAR models, with lags, to the series. Initialize each fit by specifying the first four observations. Store the Akaike information criteria (AIC) of the fits.
T = size(tbl3,1); % Total sample size
numseries = size(tbl3,2);
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
Y0 = tbl3{idxpre,:}; % Presample
Y = tbl3{idxest,:};  % Estimation sample
for j = 1:numel(numlags)
    Mdl = varm(numseries,numlags(j));
    Mdl.SeriesNames = tbl3.Properties.VariableNames;
    EstMdl(j) = estimate(Mdl,Y,'Y0',Y0);
    results = summarize(EstMdl(j));
    aic(j) = results.AIC;
end

[~,bestidx] = min(aic);
% Select the model that yields the best fit
BestMdl = EstMdl(bestidx);

resid2 = infer(BestMdl,table2array(tbl3));% calculate the residuals when using all data

Tideht_yhat2 = tbl3.TideHeight((length(tbl3.pCO2)-length(resid2)+1):end) - resid2(:,1);% estimate tide ht
Salinity_yhat2 = tbl3.Salinity((length(tbl3.pCO2)-length(resid2)+1):end) - resid2(:,2);% estimate sal
Turbidity_yhat2 = tbl3.Turbidity((length(tbl3.pCO2)-length(resid2)+1):end) - resid2(:,3);% estimate turb
Temp_yhat2 = tbl3.Temp((length(tbl3.pCO2)-length(resid2)+1):end) - resid2(:,4);% estimate temp
yhat2 = tbl3.pCO2((length(tbl3.pCO2)-length(resid2)+1):end) - resid2(:,5);% estimate pCO2

% X = table(Tideht_yhat,Salinity_yhat.*10,Turbidity_yhat.*10,Temp_yhat.*10,100*yhat);
X_all = table(Tideht_yhat2,Salinity_yhat2.*10,Turbidity_yhat2.*10,Temp_yhat2.*10,100*yhat2);
% save('VARresultstemp.mat','X','X_all')
% save('HuntDockVARM_hourlytemp.mat','numlags','nummdls','resid2','tbl3','aic','BestMdl')
co2lampcropped = 100.*tbl3.pCO2((length(tbl3.pCO2)-length(resid2)+1):end);
[R,P] = corrcoef(100*yhat2(~isnan(yhat2)),(co2lampcropped(~isnan(yhat2)))) % How well VAR fits pCO2 obs

%% Subplots
fig = figure();
colorsnew = colormap(parula(15));
eventimes = co2lamp_filteredtimes(2:2:end);

tiledlayout(4,1,'TileSpacing','none','Padding','Compact');

nexttile
% plot(4,1,1)
hold on
plot(co2lamp_filteredtimes,co2lamp_filtered,'color',colorsnew(13,:),'LineWidth',2)
plot(co2lamp_filteredtimes((length(tbl.pCO2)-length(resid2)+1):end),100*yhat2,'.k','LineWidth',2)
hold off
legend('Observed','Modeled')
ylabel('{\it p}CO_2')% (ppm)
set(gca,'FontSize',17,'YTickLabel',1000:2000:6000,'YTick',1000:2000:6000)
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
set(gca,'xticklabel',[])
xticks(hightime)

nexttile
% plot(4,1,2)
hold on
plot(data8172021.timestamp,data8172021.ECcorrected,'color',colorsnew(2,:),'LineWidth',2)
plot(co2lamp_filteredtimes,TSpCond*1000,'color',[colorsnew(3,:) 0.6],'LineWidth',2)
hold off
legend('Smart Rock','Hunt Dock')
ylabel('salinity')% (\muS/cm), is there a limit on the salinity sensor?
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
set(gca,'FontSize',17)
set(gca,'xticklabel',[],'YTickLabel',[0 10000 20000 30000 40000],'YTick',[0 10000 20000 30000 40000])
xticks(hightime)

voltage = (data8172021.turbidity2.*3.3)/4096;

right_color = colorsnew(7,:);
left_color = colorsnew(5,:);
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

nexttile
% plot(4,1,3)
yyaxis left
plot(data8172021.timestamp,voltage,'color',colorsnew(5,:),'LineWidth',2)
ylim([0 4])
ylabel('turbidity')% V
set(gca,'YColor',colorsnew(5,:),'YTickLabel',1:3,'YTick',1:3)
yyaxis right
plot(co2lamp_filteredtimes(Tturb<50),Tturb(Tturb<50),'color',[colorsnew(7,:) 0.4],'LineWidth',2)% NTU
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
xticks(hightime)
set(gca,'xticklabel',[],'YColor',[colorsnew(7,:) 0.4],'FontSize',17,'YTickLabel',[0 25 50],'YTick',[0 25 50])
legend('Smart Rock','Hunt Dock')

nexttile
% plot(4,1,4)
yyaxis left
hold on
plot(tides.DateTimeStamp,tides.cDepth,'Color',colorsnew(6,:),'LineWidth',2)
plot(co2lamp_filteredtimes,phase./25,'--','Color',[0.2627    0.1882    0.5686],'LineWidth',2)
set(gca,'YTickLabel',[0 2 4],'YTick',[0 2 4],'YColor',colorsnew(6,:))
hold off
ylabel('height')%(m)
yyaxis right
plot(co2lamp_filteredtimes,Ttemp,'s','Color',[0.0667    0.1529    0.3216 0.4])
legend('tide ht','moon phase','temp')
xlim([co2lamp_filteredtimes(1),co2lamp_filteredtimes(end)])
ylabel('temp')%(deg C)
xtickformat('D')
xticks(hightime)
xlabel('Day of Year')
set(gca,'FontSize',17,'YColor',[0.0667    0.1529    0.3216 0.4],'YTick',[27 30 33],'YTickLabel',[27 30 33])

%% Calculating net import/export of CO2
% Subtract atmospheric CO2 concentration from pCO2
co2lamp_filtered = co2lamp_filtered - CO2atm;
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
mean(net_import_Gc_time,'omitnan')% negative result implies net tidal export
sum(net_import_Gc_time,'omitnan')% total for entire study period

mean(import,'omitnan')
mean(export,'omitnan')