function obs = initObsStruct(p, num)
MAXPRNGPS = p.gps.max_prn; MAXPRNGLO = p.glo.max_prn;
MAXPRNGAL = p.gal.max_prn; MAXPRNBDS = p.bds.max_prn;
% Initialize obs struct
    obs.gps = struct;obs.gal = struct;obs.glo = struct;obs.bds = struct;
    obs.gps(1).type = 'L1 C/A';   obs.gps(2).type = 'L2';
    obs.gps(1).data.P = zeros(MAXPRNGPS,num);obs.gps(2).data.P = zeros(MAXPRNGPS,num);
    obs.gps(1).data.C = zeros(MAXPRNGPS,num);obs.gps(2).data.C = zeros(MAXPRNGPS,num);
    obs.gps(1).data.D = zeros(MAXPRNGPS,num);obs.gps(2).data.D = zeros(MAXPRNGPS,num);
    obs.gps(1).data.S = zeros(MAXPRNGPS,num);obs.gps(2).data.S = zeros(MAXPRNGPS,num);
    obs.gps(3).type = 'L1 P';     obs.gps(4).type = 'L2 P';
    obs.gps(3).data.P = zeros(MAXPRNGPS,num);obs.gps(4).data.P = zeros(MAXPRNGPS,num);
    obs.gps(3).data.C = zeros(MAXPRNGPS,num);obs.gps(4).data.C = zeros(MAXPRNGPS,num);
    obs.gps(3).data.D = zeros(MAXPRNGPS,num);obs.gps(4).data.D = zeros(MAXPRNGPS,num);
    obs.gps(3).data.S = zeros(MAXPRNGPS,num);obs.gps(4).data.S = zeros(MAXPRNGPS,num);
    obs.glo(1).type = 'L1 C/A';   obs.glo(2).type = 'L2 C/A';
    obs.glo(1).data.P = zeros(MAXPRNGLO,num);obs.glo(2).data.P = zeros(MAXPRNGLO,num);
    obs.glo(1).data.C = zeros(MAXPRNGLO,num);obs.glo(2).data.C = zeros(MAXPRNGLO,num);
    obs.glo(1).data.D = zeros(MAXPRNGLO,num);obs.glo(2).data.D = zeros(MAXPRNGLO,num);
    obs.glo(1).data.S = zeros(MAXPRNGLO,num);obs.glo(2).data.S = zeros(MAXPRNGLO,num);
    obs.glo(3).type = 'L1 P';     obs.glo(4).type = 'L2 P';
    obs.glo(3).data.P = zeros(MAXPRNGLO,num);obs.glo(4).data.P = zeros(MAXPRNGLO,num);
    obs.glo(3).data.C = zeros(MAXPRNGLO,num);obs.glo(4).data.C = zeros(MAXPRNGLO,num);
    obs.glo(3).data.D = zeros(MAXPRNGLO,num);obs.glo(4).data.D = zeros(MAXPRNGLO,num);
    obs.glo(3).data.S = zeros(MAXPRNGLO,num);obs.glo(4).data.S = zeros(MAXPRNGLO,num);
    obs.gal(1).type = 'E1';       obs.gal(2).type = 'E5b';
    obs.gal(1).data.P = zeros(MAXPRNGAL,num);obs.gal(2).data.P = zeros(MAXPRNGAL,num);
    obs.gal(1).data.C = zeros(MAXPRNGAL,num);obs.gal(2).data.C = zeros(MAXPRNGAL,num);
    obs.gal(1).data.D = zeros(MAXPRNGAL,num);obs.gal(2).data.D = zeros(MAXPRNGAL,num);
    obs.gal(1).data.S = zeros(MAXPRNGAL,num);obs.gal(2).data.S = zeros(MAXPRNGAL,num);
    obs.bds(1).type = 'B1';       obs.bds(2).type = 'B2b';
	obs.bds(1).data.P = zeros(MAXPRNBDS,num);obs.bds(2).data.P = zeros(MAXPRNBDS,num);
    obs.bds(1).data.C = zeros(MAXPRNBDS,num);obs.bds(2).data.C = zeros(MAXPRNBDS,num);
    obs.bds(1).data.D = zeros(MAXPRNBDS,num);obs.bds(2).data.D = zeros(MAXPRNBDS,num);
    obs.bds(1).data.S = zeros(MAXPRNBDS,num);obs.bds(2).data.S = zeros(MAXPRNBDS,num);
    
    obs.tr_prime = NaN(6,num); obs.tr_sow = NaN(1,num); obs.tr_week = NaN(1,num);
    obs.datetime = datetime.empty; obs.tr_posix = NaN(1,num);
end