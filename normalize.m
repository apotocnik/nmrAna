function [nY ind] = normalize(Y)

[val ind] = max(abs(Y));
nY = Y/val;