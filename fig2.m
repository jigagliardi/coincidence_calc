%% Is probability of calcite formation near climate forcing peaks better than random?
% This version exists to fix an artifact from previous versions and improve
% the clarity of the code.

 clear all
 close all

 %% Load data

load('fig_2_v9.mat'); % this version removes everything except EDC, NGRIP & calcite age data
%load('EDC2023.mat'); % has updated EDC dD from AICC2023 chronology

%% Load variables

% time range for plotting
t1 = 10;
t2 = 235;

% time range for random data: 
a = t1; %ka
b = t2; % ka

% EDC uncertainty
std2 = 1; %ka?

% y-axis location to plot ages as circles
meas(:,3) = 26;

% moving average filter window size; 
ka = 10;

% number of simulations in the monte carlo model
nsims = 10^4; 

%% Find peaks in Antarctic temperature 

%create time vectors
EDC_t2 = (t1:.1:t2); % time vector in 100 year increments
EDC_time_temp = EDCd2HAICC2012chron.AgekaBP; % raw EDC time vector 

% find indices for t1 and t2 in EDC data
[~,EDC_i] = min(abs(EDC_time_temp-t2)); % find index of time closest to t2
[~,EDC_i2] = min(abs(EDC_time_temp-t1)); % find index of time closest to t1

% created clipped versions of EDC data from t1 to t2
EDC_t = EDCd2HAICC2012chron.AgekaBP(EDC_i2:EDC_i); % raw EDC time vector clipped to desired interval 
EDC_c = EDCd2HAICC2012chron.DSMOW(EDC_i2:EDC_i); % raw EDC d2H vector clipped to desired interval
EDC_C = interp1(EDC_t,EDC_c,EDC_t2,'linear'); % interpolation of EDC d2H vector to 100 year resolution
% now EDC_t2 and EDC_C are the vectors we'll use for time and d2H respectively

% create version of EDC with orbital signal subtracted out to leave signal of millenial scale variation
wind = ka* 10; % moving average filter. wind = window size. 1 increment = 100 years. take years in ka * 10, to get increments
EDC_orbsig = smooth(EDC_C,wind) ; % EDC orbital signal

EDC_mill = zeros(length(EDC_C),1); % EDC millennial signal
for i = 1:length(EDC_C)
EDC_mill(i) = EDC_C(i)- EDC_orbsig(i);
end
% note to self: EDC_SG is now EDC_orbsig and EDC_DSG is now EDC_mill and SG_locs is now EDCm_locs

% 1) calc fractional uncertainty of EDC at each time step. EDC_unc/EDC_c
% 2) calcualte absolute uncertainty for EDC filtered using fractional unc.


% create smoothed versions of EDC and EDC_mill
EDC_sm = smooth(EDC_C,10);
EDC_mill_sm = smooth(EDC_mill,10);
r = 2;

pkprm = 9; % set value for minimum peak prominence
[~,EDC_locs] = findpeaks(EDC_C,'MinPeakProminence',pkprm); % creates a vector with the locations of peaks in EDC_C
[~,EDCm_locs] = findpeaks(EDC_mill,'MinPeakProminence',pkprm); % creates a vector with the locations of peaks in EDC_mill

%% interpolate EDC age uncertainty
% create uncertainty vector that matches length and time extent of EDC_t
EDC_time_temp2 = EDCAICC2012chron.AgekaBP; % raw EDC time vector from version with age uncertainties
[~,EDC_i3] = min(abs(EDC_time_temp2-t2)); % find index of time closest to t4
[~,EDC_i4] = min(abs(EDC_time_temp2-t1)); % find index of time closest to t3

EDC_UNC = EDCAICC2012chron.AgeStdDev(EDC_i4:EDC_i3);
EDC_age = EDCAICC2012chron.AgekaBP(EDC_i4:EDC_i3);

EDC_sig2 = interp1(EDC_age,EDC_UNC,EDC_t2);
EDC_sig2 = EDC_sig2(:)/std2;
EDC_sig2(1) = 0;
EDC_sig2(end) = 0;
EDC_sig2(end-1) = 0;

%% count forcing events, collate timing of them.

% number of EDC highs
EDC_cnt = length(EDC_locs) % count of EDC warm peaks
EDC_peaks = EDC_t2(EDC_locs); %ka
EDC_unc = EDC_sig2(EDC_locs); %ka

% number of EDC_mill highs
EDCm_cnt = length(EDCm_locs) % count of EDC millennial warm peaks
EDCm_peaks = EDC_t2(EDCm_locs); %ka
EDCm_unc = EDC_sig2(EDCm_locs); %ka


%% create vectors to store coincidence calculation results 

Co_obs_mst =   zeros(2,1); % coincidence between data and actual warm peaks
Co_dist_mst = zeros(2,nsims); % coincidence between data and actual warm peaks 
sig_mst= zeros(2,1); % standard deviation of synthetic coincidences
mn_mst = zeros(2,1); % mean of synthetic coincidences

%% define calcite ages

dates_mu =  meas(:,1); % calcite ages
dates_sigma=  meas(:,2)/2;% = age unc 1 sig.

%% Calculate coincidences

% Calculate coincidence between calcite ages and EDC warm peaks
warm_peaks_mu = EDC_t2(EDC_locs);
warm_peaks_sigma =  EDC_sig2(EDC_locs);
  coincidence = 0;
for i = 1:length(dates_mu)
for j= 1:length(warm_peaks_mu)
    sigma = sqrt(dates_sigma(i)^2 + warm_peaks_sigma(j)^2);
    coincidence = coincidence + normpdf(dates_mu(i),  warm_peaks_mu(j),sigma);
end
end
 coincidence_obs = coincidence / length(dates_mu);

% Calculate coincidence between calcite ages and synthetic warm peaks
coincidence_dist = zeros(nsims,1);
syn_size = length(EDC_locs);
for k = 1:nsims
  syn_peaks_mu = a + (b-a) .* rand(syn_size,1); %random times size of EDC peaks
    for h = 1:length(syn_peaks_mu)
        [M,I] =min(abs(syn_peaks_mu(h)- EDC_t2));
        syn_peaks_sigma = zeros(1,EDC_cnt);
        syn_peaks_sigma(h) = EDC_sig2(I);
    end
coincidence_syn = 0;
    for i = 1:length(dates_mu)
        for j= 1:length(syn_peaks_mu)
         sigma_syn = sqrt(dates_sigma(i)^2 + syn_peaks_sigma(j)^2);
         coincidence_syn = coincidence_syn + normpdf(dates_mu(i),  syn_peaks_mu(j),sigma_syn);
        end
    end
     coincidence_dist(k,1) = coincidence_syn / length(dates_mu);
end

Co_obs_mst(1) = coincidence_obs;
Co_dist_mst(1,:) = coincidence_dist(:);
sig_mst(1) = std(coincidence_dist);
mn_mst(1) = mean(coincidence_dist);

% Calculate coincidence between calcite ages and EDC millenial signal (EDC smooth subtracted)
warm_peaks_mu = EDCm_peaks;
warm_peaks_sigma = EDCm_unc;
  coincidence = 0;
for i = 1:length(dates_mu)
for j= 1:length(warm_peaks_mu)
    sigma = sqrt(dates_sigma(i)^2 + warm_peaks_sigma(j)^2);
    coincidence = coincidence + normpdf(dates_mu(i),  warm_peaks_mu(j),sigma);
end
end

 coincidence_obs = coincidence / length(dates_mu);

% Calculate coincidence with synthetic data
coincidence_dist = zeros(nsims,1);
syn_size = EDCm_cnt;
aa = max(meas(:,2));
bb = min(meas(:,2));

for k= 1:nsims
  syn_peaks_mu = a + (b-a) .* rand(syn_size,1);%

   for h = 1:length(syn_peaks_mu)
        [M,I] =min(abs(syn_peaks_mu(h)- EDC_t2));
        syn_peaks_sigma(h) = EDC_sig2(I);
    end

coincidence_syn = 0;
    for i = 1:length(dates_mu)
        for j= 1:length(syn_peaks_mu) % this may be the problem
         sigma_syn = sqrt(dates_sigma(i)^2 + syn_peaks_sigma(j)^2);
         coincidence_syn = coincidence_syn + normpdf(dates_mu(i),  syn_peaks_mu(j),sigma_syn);
        %  if isnan(coincidence_syn)
        %  pause
        %  end
         end
    end
     coincidence_dist(k,1) = coincidence_syn / length(dates_mu);

    
end

Co_obs_mst(2) =   coincidence_obs;
Co_dist_mst(2,:) = coincidence_dist(:);
sig_mst(2) = std(coincidence_dist);
mn_mst(2) = mean(coincidence_dist);


%% Calculate p-value
P = zeros(1,length(Co_obs_mst));
for i = 1:length(Co_obs_mst)
P(i) = 1 - normcdf(Co_obs_mst(i),mn_mst(i),sig_mst(i)) % not sure how this makes p-value
end

%% plot 

figure (1)

subplot(2,1,1) 
% plot EDC (raw and smoothed) with orbital signal not filtered out
plot(EDC_t2,EDC_C,Color=[0.8 0.8 0.8])
hold on
plot(EDC_t2,EDC_sm,'k');
hold on
plot(EDC_t2(EDC_locs),EDC_sm(EDC_locs),'ro')

% plot calcite ages as open circles above graph
plot(meas(:,1),(max(EDC_C)),'ko')%,'LineWidth',1.5)

% plot calcite ages and errors as blue bars
 for i = 1:length(meas(:,1))
 x1 = meas(i,1)+(meas(i,2)/1);
 x2 = meas(i,1)-(meas(i,2)/1);
 y1 = min(EDC_C);
 y2 = max(EDC_C);
 hold on
 pgon = polyshape([x1 x2 x2 x1],[y1 y1 y2 y2]);
 hold on
 plot(pgon,'Edgecolor','none','Facecolor',[0.7 0.8 0.9])
 hold on
 end

 % format axes
hold on
ylabel(['EDC \deltaD','(',char(8240),')'])
set(gca,'XMinorTick','on')
set(gca, 'xtick', t1:20:t2);
xlabel('Time (ka)')
box 'off'
xlim([t1 t2])
ylim([y1 y2])

% add histogram of coincidences
hold on
subplot(2,1,2)
histogram(Co_dist_mst(1,:))
xline(Co_obs_mst(1),'r','LineWidth',3)
xlabel('Coincidence');
hold off

%% 
figure (2)
% plot EDC
%plot(EDC_t2,EDC_C,Color=[0.8 0.8 0.8])
hold on
plot(EDC_t2,EDC_sm,'k');
hold on

% plot AIMs as red circles
plot(EDC_t2(EDC_locs),EDC_sm(EDC_locs),'ro')

% plot calcite ages and 2(?) sigma errors as blue boxes
 % for i = 1:length(meas(:,1))
 % x1 = meas(i,1)+(meas(i,2)/1);
 % x2 = meas(i,1)-(meas(i,2)/1);
 % y1 = -450;
 % y2 = -360;
 % hold on
 % pgon = polyshape([x1 x2 x2 x1],[y1 y1 y2 y2]);
 % hold on
 % plot(pgon,'Edgecolor','none','Facecolor',[0.7 0.8 0.9])
 % hold on
 % end

% format axes
box 'off'
xlim([0 120])
ylim([-450 -360])
ylabel(['EDC \deltaD','(',char(8240),')'])
xlabel('Time (ka)')
hold on
set(gca, 'xtick', t1:20:t2)
set(gca,'XMinorTick','on')
hold on