function [p, eph, obs] = readDataFiles(p,files)
%-------------------------------------------------------------------------%
p.mode_sps = 0;
p.mode_ppp = 1;
p.mode_dgnss = 2;
p.mode_rtkfloat = 3;
p.mode_rtkfix = 4;

%Define the constant parameters for GNSS system
%-------------------------------------------------------------------------%
p.gps.max_prn = 33; % The max number of GPS satellites
p.gal.max_prn = 38; % The max number of GAL satellites
p.glo.max_prn = 34; % The max number of GLO satellites
p.bds.max_prn = 65; % The max number of BDS satellites

% Get ephemeris data (.nav file, RINEX verion 3.03)
eph = parserGnssEph(p, files.eph);
% Get observables data (.obs file, RINEX verion 3.03)
obs = parserGnssObs(p, files.obs);
p.t = datetime(obs.tr_prime');

p.L2enable = false;

p.enableGPS = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO = 0; % Enable GLO: 1 means enable, 0 means close
p.enableGAL = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS = 1; % Enable BDS: 1 means enable, 0 means close
p.IGS_enable = 1; % Enable IGS correction: 1 means enable, 0 means close
p.bia_type = 0; % 0 means using CNE code bias

end