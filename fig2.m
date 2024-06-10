clear all
close all

load('fig2.mat');

%% define values

%time range for plotting
t1 = 0;
t2 = 235;

% EDC unc
std2 = 1;

%time range for random data: 
a = 0; %ka
b = 260; % ka

meas(:,3)=26;

%MA filter window size; 
ka = 10;
nsims = 10^4;
ast_unc = 0.1; %ka
%% find peaks in NH and SH insolation (precession)

NHSI_t = InsolationData5Ma.Timeka(1:401);
SHSI_t = NHSI_t;
NHSI_c = InsolationData5Ma.Insolation65NJul(1:401);
SHSI_c = InsolationData5Ma.Insolation65SJan(1:401);
NHSI_n = normalize(NHSI_c ,'zscore');
SHSI_n = normalize(SHSI_c ,'zscore');

[NHSI_pks,NHSI_locs] = findpeaks(NHSI_n,'MinPeakProminence',0.2);
[SHSI_pks,SHSI_locs] = findpeaks(SHSI_n,'MinPeakProminence',0.2);

%% find peaks in Integrated summer energy (e.g. summer length/obliquity)
ISE_t = Integ_sum_eng(1:401,1);
ISE_c = Integ_sum_eng(1:401,2);
ISE_n =normalize(Integ_sum_eng(1:401,2),'zscore');

[ISE_pks,ISE_locs] = findpeaks(ISE_n,'MinPeakProminence',0.1);

%% Find peaks in Antarctic temperature 
EDC_t = EDCd2HAICC2012chron.AgekaBP(600:4917);% 
EDC_c = EDCd2HAICC2012chron.DSMOW(600:4917);

%% moving average filter
EDC_t2 = (0:0.1:260);
EDCt = EDCd2HAICC2012chron.AgekaBP(1:4917);% 
EDCc = EDCd2HAICC2012chron.DSMOW(1:4917);
EDC_C = interp1(EDCt,EDCc,EDC_t2,'linear');

wind = ka* 10; % 1increment = 100 years. take years in ka * 10, to get increments

EDC_SG = smooth(EDC_C,wind) ; % 1integer = 100 years

%subtract out longterm trend to leave millenial scale variation

for i = 1:length(EDC_C)
EDC_DSG(i) = EDC_C(i)- EDC_SG(i);
end

%SG_locs = [SG_locs1, SG_loct];

%SG_locs = SG_locs_new;
%% EDC age unc interpolate

%create unc vector that matches length and time extent as EDC_t

EDC_unc = EDCAICC2012chron.AgeStdDev(1:1645);
EDC_age = EDCAICC2012chron.AgekaBP(1:1645);

EDC_sig = interp1(EDC_age,EDC_unc,EDC_t);
EDC_sig = EDC_sig(:)/std2;

EDC_sig2 = interp1(EDC_age,EDC_unc,EDC_t2);
EDC_sig2 = EDC_sig2(:)/std2;

%%
EDC_sm = smooth(EDC_C,10);
EDC_DSG_sm = smooth(EDC_DSG,10);
r = 2;

%% count forcing events, collate timing of them.

% number of NHSI highs
NHSI_cnt = length(NHSI_pks); % count of 
%NHSI_cnt = length(meas(:,1)); % count of 
NHSI_dates = NHSI_t(NHSI_locs); %ka
%for i =1:length(meas(:,1))
%NHSI_dates(i) = meas(i,1) + (1* rand(1)); %ka
%end
NHSI_unc = NHSI_dates; %ka
NHSI_unc(:) = ast_unc;

% number of SHSI highs
SHSI_cnt = length(SHSI_pks); % count of 
SHSI_dates = SHSI_t(SHSI_locs); %ka
SHSI_unc = SHSI_t(SHSI_locs); %ka
SHSI_unc(:) = ast_unc;

% number of ISE highs
ISE_cnt = length(ISE_pks); % count of 
ISE_dates = ISE_t(ISE_locs); %ka
ISE_unc = (ISE_t(ISE_locs)); %ka
ISE_unc(:) = ast_unc;

% number of EDC highs
EDC_cnt = length(SG_locs); % count of
EDC_dates = EDC_t2(SG_locs); %ka
EDC_unc = EDC_sig2(SG_locs); %ka

% number of EDC-LR04 highs
SG_cnt = length(SG_locs); % count of 
SG_dates = EDC_t2(SG_locs); %ka
SG_unc = EDC_sig2(SG_locs); %ka

cnt_master(1) = NHSI_cnt;
cnt_master(2) = SHSI_cnt;
cnt_master(3) = ISE_cnt;
cnt_master(4) = EDC_cnt;
cnt_master(5) = SG_cnt;

Co_obs_mst =   zeros(5,1);
Co_dist_mst = zeros(5,nsims);
sig_mst= zeros(5,1);
mn_mst = zeros(5,1);
%% compare the forcing peaks with the MC synthetic ages
dates_mu =  meas(:,1); % calcite ages
dates_sigma=  meas(:,2)/2;% = age unc 1 sig.

%% Calculate coincidences

 % the product of two normal PDFs is a curve with an integrated area equal to another normal PDF,
 % specifically, ∫ pdf(Normal(μ₁, σ₁), x) * pdf(Normal(μ₂, σ₂), x) from x = -inf to inf
 % = normpdf(μ₁, sqrt(σ₁^2 + σ₂^2), μ₂)
 %This lets us calculate coincidence products without numerical integration

% NHSI 
warm_peaks_mu = NHSI_dates;
warm_peaks_sigma = NHSI_unc;
  coincidence = 0;
for i = 1:length(dates_mu)
for j= 1:length(warm_peaks_mu)
    sigma = sqrt(dates_sigma(i)^2 + warm_peaks_sigma(j)^2);
    coincidence = coincidence + normpdf(dates_mu(i), warm_peaks_mu(j),sigma);
end
end

 coincidence_obs = coincidence / length(dates_mu);

% Calculate coincidence with synthetic data
coincidence_dist = zeros(nsims,1);
syn_size = NHSI_cnt;
syn_peaks_sigma = NHSI_unc;

for k= 1:nsims
  syn_peaks_mu = a + (b-a) .* rand(syn_size,1);%random times size of NHSI peaks

coincidence_syn = 0;
    for i = 1:length(dates_mu)
        for j= 1:length(syn_peaks_mu)
         sigma_syn = sqrt(dates_sigma(i)^2 + syn_peaks_sigma(j)^2);
         coincidence_syn = coincidence_syn + normpdf(dates_mu(i),syn_peaks_mu(j),sigma_syn);
        end
    end
     coincidence_dist(k,1) = coincidence_syn / length(dates_mu);
end

Co_obs_mst(1) =   coincidence_obs;
Co_dist_mst(1,:) = coincidence_dist(:);
sig_mst(1) = std(coincidence_dist);
mn_mst(1) = mean(coincidence_dist);

% SHSI 
warm_peaks_mu = SHSI_dates;
warm_peaks_sigma = SHSI_unc;
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
syn_size = SHSI_cnt;
syn_peaks_sigma = SHSI_unc;

for k= 1:nsims
  syn_peaks_mu = a + (b-a) .* rand(syn_size,1);%random times size of NHSI peaks
% 
coincidence_syn = 0;
    for i = 1:length(dates_mu)
        for j= 1:length(syn_peaks_mu)
         sigma_syn = sqrt(dates_sigma(i)^2 + syn_peaks_sigma(j)^2);
         coincidence_syn = coincidence_syn + normpdf(dates_mu(i),  syn_peaks_mu(j),sigma_syn);
        end
    end
     coincidence_dist(k,1) = coincidence_syn / length(dates_mu);
end

Co_obs_mst(2) =   coincidence_obs;
Co_dist_mst(2,:) = coincidence_dist(:);
sig_mst(2) = std(coincidence_dist);
mn_mst(2) = mean(coincidence_dist);

% ISE
warm_peaks_mu = ISE_dates;
warm_peaks_sigma = ISE_unc;
  coincidence = 0;
for i = 1:length(dates_mu)
for j= 1:length(warm_peaks_mu)
    sigma = sqrt(dates_sigma(i)^2 + warm_peaks_sigma(j)^2);
    coincidence = coincidence + normpdf(dates_mu(i), warm_peaks_mu(j), sigma);
end
end

 coincidence_obs = coincidence / length(dates_mu);

% Calculate coincidence with synthetic data
coincidence_dist = zeros(nsims,1);
syn_size = ISE_cnt;
syn_peaks_sigma = ISE_unc;

for k= 1:nsims
  syn_peaks_mu = a + (b-a) .* rand(syn_size,1);%random times size of NHSI peaks


coincidence_syn = 0;
    for i = 1:length(dates_mu)
        for j= 1:length(syn_peaks_mu)
         sigma_syn = sqrt(dates_sigma(i)^2 + syn_peaks_sigma(j)^2);
         coincidence_syn = coincidence_syn + normpdf(dates_mu(i),  syn_peaks_mu(j),sigma_syn);
        end
    end
     coincidence_dist(k,1) = coincidence_syn / length(dates_mu);
end

Co_obs_mst(3) =   coincidence_obs;
Co_dist_mst(3,:) = coincidence_dist(:);
sig_mst(3) = std(coincidence_dist);
mn_mst(3) = mean(coincidence_dist);

% EDC 
warm_peaks_mu = EDC_t2(SG_locs);
warm_peaks_sigma =  EDC_sig2(SG_locs);
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
syn_size = length(SG_locs);

for k= 1:nsims
  syn_peaks_mu = a + (b-a) .* rand(syn_size,1);%random times size of NHSI peaks

    for h = 1:length(syn_peaks_mu)
        [M,I] =min(abs(syn_peaks_mu(h)- EDC_t2));
        syn_peaks_sigma(h) = EDC_sig2(I);
    end
% 
coincidence_syn = 0;
    for i = 1:length(dates_mu)
        for j= 1:length(syn_peaks_mu)
         sigma_syn = sqrt(dates_sigma(i)^2 + syn_peaks_sigma(j)^2);
         coincidence_syn = coincidence_syn + normpdf(dates_mu(i),  syn_peaks_mu(j),sigma_syn);
        end
    end
     coincidence_dist(k,1) = coincidence_syn / length(dates_mu);
end

Co_obs_mst(4) =   coincidence_obs;
Co_dist_mst(4,:) = coincidence_dist(:);
sig_mst(4) = std(coincidence_dist);
mn_mst(4) = mean(coincidence_dist);

% EDC Millenial  
warm_peaks_mu = SG_dates;
warm_peaks_sigma = SG_unc;
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
syn_size = SG_cnt;
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
        for j= 1:length(syn_peaks_mu)
         sigma_syn = sqrt(dates_sigma(i)^2 + syn_peaks_sigma(j)^2);
         coincidence_syn = coincidence_syn + normpdf(dates_mu(i),  syn_peaks_mu(j),sigma_syn);
        end
    end
     coincidence_dist(k,1) = coincidence_syn / length(dates_mu);
end

Co_obs_mst(5) =   coincidence_obs;
Co_dist_mst(5,:) = coincidence_dist(:);
sig_mst(5) = std(coincidence_dist);
mn_mst(5) = mean(coincidence_dist);


%% P value calc
for i = 1:length(Co_obs_mst)
P(i) =1-normcdf(Co_obs_mst(i),mn_mst(i),sig_mst(i))
end

%% plot 
figure (1)
% subplot(3,1,1)
% plot(EDC_t2,EDC_C,Color=[0.8 0.8 0.8])
% hold on
% plot(EDC_t2,EDC_sm,'k');
% hold on
% %plot(EDC_t2(SG_locs),EDC_sm(SG_locs),'ro')
% ylabel(['EDC \deltaD','(',char(8240),')'])
% box 'off'
% xlim([t1 t2])
% 
% set(gca,'xticklabel',{[]})
% set(gca, 'xtick', t1:20:t2);
% set(gca,'XMinorTick','on')
% hold off

subplot(2,1,1)
% orbital signal filtered out
%plot(EDC_t2,EDC_DSG,Color=[0.8 0.8 0.8]);
%hold on
%plot(EDC_t2,EDC_DSG_sm,'k');
%hold on
%plot(EDC_t2(SG_locs),EDC_DSG_sm(SG_locs),'ro')

% orbital signal not filtered out
plot(EDC_t2,EDC_C,Color=[0.8 0.8 0.8])
hold on
plot(EDC_t2,EDC_sm,'k');
hold on
plot(EDC_t2(SG_locs),EDC_sm(SG_locs),'ro')

% plot calcite ages as open circles above graph
plot(meas(:,1),-350,'ko')%,'LineWidth',1.5)

 for i = 1:length(meas(:,1))
 x1 = meas(i,1)+(meas(i,2)/1);
 x2 = meas(i,1)-(meas(i,2)/1);
 y1 = -450;
 y2 = -360;
 hold on
 pgon = polyshape([x1 x2 x2 x1],[y1 y1 y2 y2]);
 hold on
 plot(pgon,'Edgecolor','none','Facecolor',[0.7 0.8 0.9])
 hold on
 end

hold on
ylabel(['EDC \deltaD','(',char(8240),')'])
set(gca,'XMinorTick','on')
set(gca, 'xtick', t1:20:t2);
xlabel('Time (ka)')
box 'off'
xlim([t1 t2])
ylim([y1 -345])

hold all
subplot(2,1,2)
%title('forcing: millenial EDC')
histogram(Co_dist_mst(5,:))
xline(Co_obs_mst(5),'r','LineWidth',3)
%yline(mean(Co_dist_mst(5,:)),'k')
xlabel('Coincidence');
% ylabel('Frequency');
%set(gca,'xticklabel',{[]})
%set(gca,'YAxisLocation','right')
%set(gca, 'YScale', 'log')
 hold off

% figure (3)
% 
% subplot(2,1,1)
% plot(1:5, P,'ko','MarkerSize',8,'markerfacecolor', 'b')
% yline(0.05,'r')
% %set(gca, 'YScale', 'log')
% hold on
% ylabel('Signifigance P-value');
% 
% subplot(2,1,2)
% plot(1:5, cnt_master,'ko','MarkerSize',8,'markerfacecolor', 'b')
% hold on
% xlabel('Experiment number');
% ylabel('Peak count');
%% 
figure (2)
% plot EDC
plot(EDC_t2,EDC_C,Color=[0.8 0.8 0.8])
hold on
plot(EDC_t2,EDC_sm,'k');
hold on

% plot AIMs as red circles
plot(EDC_t2(SG_locs),EDC_sm(SG_locs),'ro')

% plot calcite ages and 2(?) sigma errors as blue boxes
 for i = 1:length(meas(:,1))
 x1 = meas(i,1)+(meas(i,2)/1);
 x2 = meas(i,1)-(meas(i,2)/1);
 y1 = -450;
 y2 = -360;
 hold on
 pgon = polyshape([x1 x2 x2 x1],[y1 y1 y2 y2]);
 hold on
 plot(pgon,'Edgecolor','none','Facecolor',[0.7 0.8 0.9])
 hold on
 end
 
% format axes
box 'off'
xlim([t1 t2])
ylim([-450 -360])
ylabel(['EDC \deltaD','(',char(8240),')'])
xlabel(['Time (ka)'])
hold on
set(gca, 'xtick', t1:20:t2)
set(gca,'XMinorTick','on')
hold on
