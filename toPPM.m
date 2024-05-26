function [ppm] = toPPM(kHz,reverse)
% If reverse == 1, then we have reverse!!! ppm -> kHz
FR = evalin('base','FR');
REF = evalin('base','REF');

if nargin < 2
    reverse = 0;
end
    
if reverse == 0
    ppm = (kHz+1000*(FR(1)-REF))/REF*1000;
else
%     kHz = ppm/1000*REF-1000*(FR(1)-REF);
%   ppm<-kHz, kHz<-ppm
    ppm = kHz/1000*REF-1000*(FR(1)-REF);
end

% if FR(1)<100 && FR(1)>90
%     REF = 95.557;   % 380 MHz 13C
% elseif FR(1)>120 && FR(1)<135
%     REF = 130.8712; % 400 MHz 87Rb
% end