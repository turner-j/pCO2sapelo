%% MVGC statistics demo
%
% Demonstrates MVGC toolbox time series statistical and spectral analysis tools.
% PLEASE NOTE: IN ORDER TO RUN THIS SCRIPT YOU MUST DOWNLOAD THE MVGC
% TOOLBOX AVAILABLE THROUGH MATLAB ONLINE AND HAVE THIS SCRIPT OPEN IN THAT
% FOLDER.
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters
clear; clc;
load('VARresultstemp.mat')
load('HuntDockVARM_hourlytemp.mat')

ntrials   = 1;      % number of trials
nobs      = height(X_all);   % number of observations

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default [set in startup.m])
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default [set in startup.m])

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 5;     % maximum model order for model order estimation

acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC: 'F' for Granger's F-test, 'chi2' for Geweke's chi2 test or leave empty for default
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 22/86400;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)
specm     = [];     % power spectrum estimation method: 'WELCH' (Welch method - default) or 'MT' (multi-taper)

etests    = false;  % do experimental (unit-root stationarity) tests?
stlags    = [];     % number of lags for stationarity tests (or leave empty for automatic default)

acorr     = true;   % plot autocorrelation (else autocovariance)?

%% Generate VAR data
% Get VAR coefficients

for i = 1:length(BestMdl.AR)
    AT(:,:,i) = cell2mat(BestMdl.AR(i));% Autoregressive coefficient matrices associated with the lagged responses
end

nvars = size(AT,1); % number of variables

% Residuals covariance matrix.
SIGT = BestMdl.Covariance;

% Call VAR time series data

X = table2array(tbl3);X = X(8:end,:);X = X.';

[row, col] = find(isnan(X));
X(row,col) = X(row,col+1);

% Inputs were originally scaled, so scale back to original
tbl3.Salinity = tbl3.Salinity.*10;
tbl3.Turbidity = tbl3.Turbidity.*10;
tbl3.Temp = tbl3.Temp.*10;
tbl3.pCO2 = tbl3.pCO2.*100;
%% Model order estimation

% Calculate information criteria up to max model order

[~,bmo_AIC] = min(aic);

% Plot information criteria.

figure(1); clf;
plot(aic);
legend('AIC');

amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
fprintf('actual model order     = %d\n',amo);

% Select model order

if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = bmo_AIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = bmo_BIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% VAR model estimation and autocovariance calculation

% Calculate VAR model; return residuals E too, since we need them later for
% statistical routines.

A = AT;% VAR coefficients matrix
SIG = SIGT;% residuals covariance matrix
E = resid2';
%E = E(:,6:end);% residuals time series

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

%% Now calculate autocovariance according to the VAR model, to as many lags
% as it takes to decay to below the numerical tolerance level, or to acmaxlags
% lags if specified (i.e. non-empty).

[G,info] = var_to_autocov(AT,SIGT,[]);

F = autocov_to_pwcgc(G);
cd = mean(F(~isnan(F)));
fprintf('\ncausal density = %f\n',cd);


%% Check for failed spectral GC calculation
f = autocov_to_spwcgc(G,fres);

assert(~isbad(f,false),'spectral GC calculation failed');

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

figure(4); clf;
% sgtitlex('Pairwise-conditional Granger causality - time domain');
subplot(1,3,1);
plot_pw(F,flipud(gray));
title('Pairwise-conditional GC');
set(gca,'FontSize',16)

subplot(1,3,2);
plot_pw(pval,[]);
title('p-values');
set(gca,'FontSize',16)

subplot(1,3,3);
plot_pw(sig,flipud(gray));
title(['Significant at p = ' num2str(alpha)])
set(gca,'FontSize',16)

% Report and check for errors.
var_acinfo(info,true); % report results (and bail out on error)

%% Spectral analysis

[S,fres] = autocov_to_cpsd(G,fres); % for model
SE = tsdata_to_cpsd(X,fres,specm);  % from data (empirical)

% plot (auto-)spectra

figure(6); clf;
plot_cpsd3(cat(4,S,SE),{'model','data'},fs,[],true);
%% VAR stats tests
[A,SIG,E] = tsdata_to_var(X,morder,regmode);
% Check that residuals are white (Durbin-Watson test).

[dw,dwpval] = whiteness(X,E);
fprintf('\nDurbin-Watson statistics =\n'); disp(dw);
dwsig = significance(dwpval,alpha,mhtc); % significance adjusted for multiple hypotheses
notwhite = find(dwsig);
if isempty(notwhite)
    fprintf('all residuals are white by Durbin-Watson test at significance %g\n',alpha);
else
    fprintf(2,'WARNING: autocorrelated residuals at significance %g for variable(s): %s\n',alpha,num2str(notwhite));
end

% Check R^2 stats.

[~,RSQADJ] = rsquared(X,E);
fprintf('\nRSQ (adjusted) =\n'); disp(RSQADJ);
rsqthreshold = 0.3; % like GCCA
badqsq = find(RSQADJ < rsqthreshold);
if isempty(badqsq)
    fprintf('adjusted r-squares OK: > %g%% of variance is accounted for by the model\n',100*rsqthreshold);
else
    fprintf(2,'WARNING: low adjusted r-square values (< %g) for variable(s): %s\n',rsqthreshold,num2str(badqsq));
end

% Check model consistency (ie. proportion of correlation structure of the data
% accounted for by the model).

cons = 100*consistency(X,E); % percent
fprintf('\nmodel consistency = %.0f%%\n',cons);
consthreshold = 80;          % like GCCA
if cons > consthreshold
    fprintf('consistency OK: > %g%%\n',consthreshold);
else
    fprintf(2,'WARNING: low consistency (< %g%%)\n',consthreshold);
end

