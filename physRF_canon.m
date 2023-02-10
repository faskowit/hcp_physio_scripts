function [HR_conv,RVT_conv,xPhys] = physRF_canon(HR,resp,samprate,BOLDind,GS)

if nargin < 3
    samprate = 10 ;
end

if nargin < 5
    GS = [] ; 
end

%%

HR = HR(:) ;
resp = resp(:) ;

Ts_srate = 1/samprate ;  % Sampling period in seconds
time_srate = 0:Ts_srate:(length(HR)-1)*Ts_srate;

[pks,loc] = findpeaks(resp,time_srate,'MinPeakDistance',2,'MinPeakHeight',0.2);
respUpp = interp1([0,loc,time_srate(end)],[pks(1),pks',pks(end)],time_srate);

[pks,loc] = findpeaks(-resp,time_srate,'MinPeakDistance',2,'MinPeakHeight',0.2); 
pks=-pks;
respLow = interp1([0,loc,time_srate(end)],[pks(1),pks',pks(end)],time_srate);

BR = 60./diff(loc);
time_BR = [time_srate(1),(loc(2:end)-loc(1:end-1))/2+loc(1:end-1),time_srate(end)];
BR = interp1(time_BR,[BR(1),BR,BR(end)],time_srate);

RVT = ((respUpp-respLow).*BR)';

t_IR = 0:Ts_srate:60;
RRF = 0.6*t_IR.^(2.1).*exp(-t_IR/1.6)-0.0023*t_IR.^3.54.*exp(-t_IR/4.25);
RRF  = RRF/max(RRF);
CRF = 0.6*t_IR.^2.7.*exp(-t_IR/1.6)-(16/(sqrt(2*pi*9)))*exp(-(t_IR-12).^2/18);
CRF  = CRF/max(CRF);

% convolve
HR_conv = pad_conv(smooth(HR,6*samprate),CRF,round(length(HR).*.3));
RVT_conv = pad_conv(RVT,RRF,round(length(RVT).*.3));

% %% chunk downsample
% % doesn't really make a difference
% 
% chunks = [BOLDind [ BOLDind(2:end)-1 ; BOLDind(end)+(BOLDind(end)-BOLDind(end-1)) ] ] ;
% HR_conv2 = zeros(size(chunks,1),1) ;
% for idx = 1:size(chunks,1)
%     ii = chunks(idx,:) ; 
%     HR_conv2(idx) = mean(HR_conv(ii(1):ii(2))) ; 
% end

%% regressors

xPhys = [HR_conv(BOLDind),RVT_conv(BOLDind)];   
xPhys = detrend(xPhys,'linear');
xPhys = zscore(xPhys) ;

%% regress against GS

if ~isempty(GS)

    NV = length(GS);
    
    regr = [ones(NV,1) xPhys];
    
    B = regr\GS;     yPred = regr*B;
    
    r_PRF(1) = corr(yPred,GS);
    yPred_card = regr(:,2)*B(2);  r_PRF(2) = corr(yPred_card,GS);
    yPred_resp = regr(:,3)*B(3);  r_PRF(3) = corr(yPred_resp,GS);
    
    fprintf(' ----------------------------------------------- \n')
    fprintf('Correlation b/w GS and PRF output \n')
    fprintf('CRF (HR): %3.2f  \n',r_PRF(2))
    fprintf('RRF (RVT): %3.2f  \n',r_PRF(3))
    fprintf('CRF & RRF (HR & RVT): %3.2f  \n',r_PRF(1))

end