function [ aa ] = physRF_retroicorr(resp,ppg,trigger,samprate)
% get retroicorr regressors, in table form

NR = 2 ; 

N = length(resp);
Fs = samprate; 
Ts=1/Fs;  
timePhys=0:Ts:(N-1)*Ts;

ret_resp = func_RETR_Resp_regressors(resp,NR,Fs) ;
ret_card = func_RETR_Card_regressors(timePhys,ppg,NR) ;

trig_diff = zeros(length(trigger),1) ;
trig_diff(diff(trigger)==1) = 1 ;
trig_diff(1) = 1 ;

ret_resp_dwnsmp = ret_resp(~~trig_diff,:) ;
ret_card_dwnsmp = ret_card(~~trig_diff,:) ;
 
aa = [ret_resp_dwnsmp ret_card_dwnsmp] ;

% make a table
% resp_names = append(repmat({'resp_ricor'},4,1),cellstr(string(1:NR*2)')) ;
% card_names = append(repmat({'card_ricor'},4,1),cellstr(string(1:NR*2)')) ;
% tab = array2table(aa) ;
% tab.Properties.VariableNames = [ resp_names ;  card_names ] ;
