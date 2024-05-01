function spec_cov_ecef = specInNedToEcef(pos_lla, spec)

% pos_lla: position in Lat, Lng, H (rad, rad, m) (3 x 1)
% spec: accuracy specification in ENU (3 x 3)
R_e2g=ll2R(pos_lla);
spec_cov_ecef = R_e2g * spec * R_e2g';