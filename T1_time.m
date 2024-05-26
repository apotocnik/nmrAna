function [cas] = T1_time(NS,TAUMIN,TAUMAX)
D0 = 0.05;
NT1 = 36;
KT1 =(TAUMAX/TAUMIN)^(1/(NT1-1));
cas = round(100*NS*(NT1*D0+TAUMIN*(1-KT1^NT1)/(1-KT1))/3600)/100;
end