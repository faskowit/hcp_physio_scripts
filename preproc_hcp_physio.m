function [HR_resmp,resp_resmp,time_resmp,trig_resmp,...
    PPGlocs,cardiac_resmp] = preproc_hcp_physio(inFile,Fs_resmp)

if nargin < 2
    Fs_resmp = 10 ;
end

% loadPhys = load('./mydata/100206_3T_rfMRI_REST1_LR_Physio_log.txt') ;
loadPhys = load(inFile) ;

resp_cln = zscore(loadPhys(:,2));
raw_cardiac = zscore(loadPhys(:,3));

% Create timeline at the frequency sampling Fs = 400 Hz that the recordings
% were acquried

N = length(resp_cln);
Fs = 400; % given by HCP
Ts=1/Fs;  
timePhys=0:Ts:(N-1)*Ts;

% setup downsample
Ts_resmp = 1/Fs_resmp; 
time_resmp = timePhys(1):Ts_resmp:timePhys(end);

trig_orig = loadPhys(:,1);
tt_dwnsmp = interp1(timePhys,trig_orig,time_resmp,"nearest");
diff_dwnsmp = diff(tt_dwnsmp)==1 ;
trig_resmp = zeros(length(tt_dwnsmp),1) ;
trig_resmp(diff_dwnsmp==1) = 1 ; 
trig_resmp(1) = 1 ; % add the first trigger

%    =========================================
%% 3: Preprocess the cardiac signal and extract the heart rate (HR)

% The cardiac signal is first band-pass filtered and, subsequently, the
% peaks in the signal are identified. Based on the times of the peaks, the
% HR is estimated. As the cardiac signal is often noisy, there are some
% steps that we can take in order to improve the quality of the extracted
% HR.

% First, we have to specify the 'minPeak' variable which denotes the
% minimum time interval in seconds between peaks needed for the identification of the
% peaks. Normal values for 'minPeak' vary from 0.50 to 0.90 seconds. This value should be
% adjusted by the user for each scan separately based on visual inspection

% Then, we have to specify some parameters related to the correction for
% outliers in the extracted HR based on visual inspection. Typically, the
% 'filloutliers_ThresholdFactor' has to be somewhere between 3-20 in order
% to correct for outliers in noisy epochs while also retain real abrupt changes
% in HR. By zooming in a time interval with an abrupt change in HR, based
% on the quality of the cardiac signal and the time of peaks we can be more
% confident whether this abrupt change is real or artefactual.

%  Set the following parameters !!

minPeak = 0.55;
filloutliers_window = 30*Fs;     % given in number of samples
filloutliers_ThresholdFactor = 15;     %  normal range between 3-20
%  ---------------------------------------

% filter signal
f1 = 0.3; f2 = 10; [filt_b,filt_a] = butter(2,[f1,f2]/(Fs/2));
cardiac_filt = filtfilt(filt_b,filt_a,raw_cardiac);
[~,PPGlocs] = findpeaks(cardiac_filt,timePhys,'MinPeakDistance',minPeak);

HR = 60./diff(PPGlocs);
time_HR = [timePhys(1),(PPGlocs(2:end)-PPGlocs(1:end-1))/2+PPGlocs(1:end-1),timePhys(end)];

% downsample
HR_raw = interp1(time_HR,[HR(1),HR,HR(end)],timePhys);
HR_filloutl = filloutliers(HR_raw,'linear','movmedian',filloutliers_window,'ThresholdFactor',filloutliers_ThresholdFactor);

HR_resmp = interp1(timePhys,HR_filloutl,time_resmp); 
HR_resmp=HR_resmp(:);

cardiac_resmp = interp1(timePhys,cardiac_filt,time_resmp); 
cardiac_resmp = cardiac_resmp(:) ;

%    =========================================
%% 4: Preprocess the respiratory signal

% The respiratory signal is first low-pass filtered and, subsequently,
% linearly detrended and corrected for outliers. The 
% parameters for the outlier correction have to be adjusted for each scan
% separately based on visual inspection.

%  Set the following parameters !!

filloutliers_window = 0.3*Fs;     % given in number of samples
filloutliers_ThresholdFactor = 0.2;     
%  ---------------------------------------

f1=0.01; f2=5; [filt_b,filt_a] = butter(2,[f1,f2]/(Fs/2));

resp_cln = detrend(resp_cln,'linear');
resptmp = filloutliers(resp_cln,'linear','movmedian',filloutliers_window,'ThresholdFactor',filloutliers_ThresholdFactor);
resp_cln = zscore(filter(filt_b,filt_a,resptmp));

resp_resmp = interp1(timePhys,resp_cln,time_resmp);
resp_resmp = resp_resmp(:) ;

% make it a col
time_resmp = time_resmp(:) ;

end

