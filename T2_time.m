function [cas] = T2_time(D0,NS)
TAUMIN = 0.000060;
TAUMAX = 0.01156;
NT2 = 24;
KT2 =(TAUMAX/TAUMIN)^(1/(NT2-1));
cas = round(100*NS*(NT2*D0+TAUMIN*(1-KT2^NT2)/(1-KT2))/3600)/100;
end
