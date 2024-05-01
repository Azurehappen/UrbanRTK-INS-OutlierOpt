function corr = beidouCodeCorr(prn, elevation)

% Elevation-dependent correction values for Beidou code measurements
% Beidou Freq. band
% B1: F2, for example C2I
% B2: F7, for example C7I
% B3: F6, for example C6I
igso_prn = [6,7,8,9,10,13,16,31,38:40,56];
meo_prn = [11,12,14,19:30,32:37,41:46,57,58];

igso.b1 = [-0.55, -0.40, -0.34, -0.23, -0.15, -0.04, 0.09, 0.19, 0.27, 0.35];
igso.b2 = [-0.71, -0.36, -0.33, -0.19, -0.14, -0.03, 0.08, 0.17, 0.24, 0.33];
igso.b3 = [-0.27, -0.23, -0.21, -0.15, -0.11, -0.04, 0.05, 0.14, 0.19, 0.32];

meo.b1 = [-0.47, -0.38, -0.32, -0.23, -0.11, 0.06, 0.34, 0.69, 0.97, 1.05];
meo.b2 = [-0.40, -0.31, -0.26, -0.18, -0.06, 0.09, 0.28, 0.48, 0.64, 0.69];
meo.b3 = [-0.22, -0.15, -0.13, -0.10, -0.04, 0.05, 0.14, 0.27, 0.36, 0.47];

tElev = floor((elevation/10.0)+0.5)+1;

% The current application only considers Beidou B1 single freq.
if ismember(prn, igso_prn)
    corr = igso.b1(tElev);
elseif ismember(prn, meo_prn)
    corr = meo.b1(tElev);
else
    warning('No match for prn: %d', prn);
end

end