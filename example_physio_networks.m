clc
clearvars

%% RRF basis

RRF_p = @(t) ( 0.6.*(t.^2.1).*exp(-t./1.6) ) - ( 0.0023.*(t.^3.54).*exp(-t./4.25) )   ;
RRF_t1 = @(t) ( -0.79.*(t.^2.1).*exp(-t./1.6) ) + ( 2.66.*(t.^1.1).*exp(-t./1.6) ) ;
RRF_t2 = @(t) ( -0.069.*(t.^2.54).*exp(-t./4.25) ) + ( 0.0046.*(t.^3.54).*exp(-t./4.25) ) ;
RRF_d1 = @(t) ( 0.16.*(t.^3.1).*exp(-t./1.6) ) ;
RRF_d2 = @(t) ( 0.00014.*(t.^4.54).*exp(-t./4.25) ) ;

RRFfuncs = { RRF_p RRF_t1 RRF_t2 RRF_d1 RRF_d2 } ;

%% CRF basis

CRF_p = @(t) ( 0.3.*(t.^2.7).*exp(-t./1.6) ) - ( 1.05.*exp( (-(t-12).^2)./18 ) ) ;
CRF_t1 = @(t) ( 1.94.*(t.^1.7).*exp(-t./1.6) ) - ( 0.45.*(t.^2.7).*exp(-t./1.6) ) ;
CRF_t2 = @(t) 0.55 .* (t-12) .* exp( (-(t-12).^2)./18 ) ;
CRF_d1 = @(t) 0.056.*(t.^3.7).*exp(-t./1.6) ;
CRF_d2 = @(t) 0.15 .* ((t-12).^2) .* exp( (-(t-12).^2)./18 ) ;

CRFfuncs = { CRF_p CRF_t1 CRF_t2 CRF_d1 CRF_d2 } ;
 
%% viz it

TT = 0:0.1:50 ;

figure
for idx = 1:5

    hold on 
    plot(TT,RRFfuncs{idx}(TT))

end
hold off

figure
for idx = 1:5

    hold on 
    plot(TT,CRFfuncs{idx}(TT))

end
hold off

%% get the phyiso

ff = '/Users/faskowitzji/joshstuff/sandbox/PRF_estimation/mydata/100307_3T_rfMRI_REST1_LR_Physio_log.txt' ;
rr = 100 ;
% GS = detrend(load('/Users/faskowitzji/joshstuff/sandbox/PRF_estimation/mydata/rfMRI_REST1_LR_GLOBAL.txt'),'linear') ;

[~,resp,tt,dd,ppg,pulse] = preproc_hcp_physio(ff,rr) ;

[a,b] = physRF_pnets(resp,pulse,rr,0.72,1200) ; 

addpath('./ext/')
tab = physRF_retroicorr(resp,ppg,dd,rr) ;

%% smooth out HR and and resp

% add some minpreakdistances?

% 0.2 seconds = 
c_min_dist = rr*0.2 / 2 ;
[c_peaks,c_loc] = findpeaks(pulse,'MinPeakProminence',prctile(abs(pulse),95).*0.5,...
    'MinPeakDistance',c_min_dist) ; 
c_peak_ind = zeros(length(pulse),1) ; 
c_peak_ind(c_loc) = 1 ;

interbeatint = get_contact_times(~c_peak_ind) ;
c_peak_intervals = nan(length(pulse),1) ;
% each beat gets a inter-beat-interval assigned to it, in seconds
c_peak_intervals(c_loc) = interbeatint(1:(end-1)) ./ rr  ;

% 0.1 seconds 
r_min_dist = rr*0.1 / 2 ;
[r_up_peaks,r_up_loc] = findpeaks(resp,...
    'MinPeakProminence',prctile(resp,95).*0.01,'MinPeakDistance',r_min_dist) ; 
[r_dn_peaks,r_dn_loc] = findpeaks(-resp,...
    'MinPeakProminence',prctile(-resp,95).*0.01,'MinPeakDistance',r_min_dist) ; 
r_dn_peaks = -r_dn_peaks ;

r_peak_amps = nan(length(resp),1) ;
r_peak_amps(r_up_loc) = r_up_peaks ;
r_peak_amps(r_dn_loc) = r_dn_peaks ;

%% make 6 second sliding window

windowsz = 6 * rr ;
ntp = length(dd)-1 ;
TR = 0.72 ;
slidesz = TR*rr ; 

% calculate how many tvMats we will calculate
winInds = (0:slidesz:ntp-1)' + (1:windowsz) ;
% adjust to center the window, lazy way?
winInds = round(winInds - floor(windowsz/2)) ; 
% calculate rows that exceed the ntp, or are less than 0
winInds(winInds>ntp) = 0 ;
winInds(winInds<1) = 0 ;

%% get signals

% cardiac mean within 6s windows
HBI_sig = arrayfun(@(ii_) mean(c_peak_intervals(nonzeros(winInds(ii_,:))),'omitnan'),1:size(winInds,1) ) ; 


% respiratory std within 6s windows
RV_sig = arrayfun(@(ii_) std(r_peak_amps(nonzeros(winInds(ii_,:))),'omitnan'),1:size(winInds,1) ) ; 

%% make the regressors 

nTR = 1200 ;

r1 = cell2mat(arrayfun(@(ii_) pad_conv(HBI_sig,RRFfuncs{ii_}(1:50)) ,1:5,'UniformOutput',false)) ;
r1 = zscore(r1(1:nTR,:)) ;

r2 = cell2mat(arrayfun(@(ii_) pad_conv(RV_sig,CRFfuncs{ii_}(1:50)) ,1:5,'UniformOutput',false)) ;
r2 = zscore(r2(1:nTR,:)) ;





