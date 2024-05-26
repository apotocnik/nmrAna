%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phaseSmooth
%
% Smooth phase vector to avoid discontinuities due to modulus
%
% Anton Potoènik, 19.12.2013, F5 @ IJS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% phi = [0 2pi)
%

function phi = phaseSmooth(phi)
borders = {1};
b = 1;

% Find borders
for i = 1:numel(phi)-1
    if abs(phi(i)-phi(i+1)) > pi, 
        b = b + 1; 
        borders{b} = i+1;
    else
        borders{b} = [borders{b} i+1];
    end
end

% Find the largest domain
[m mind] = max(cellfun(@numel,borders));

% smooth to the left
for bi=mind-1:-1:1
    dif = (phi(borders{bi}(end)) - phi(borders{bi+1}(1)))/2/pi;
    n = round(dif);
    phi(borders{bi}) = phi(borders{bi}) - n*2*pi;
end

% smooth to the right
for bi=mind+1:b
    dif = (phi(borders{bi-1}(end)) - phi(borders{bi}(1)))/2/pi;
    n = round(dif);
    phi(borders{bi}) = phi(borders{bi}) + n*2*pi;
end

