function eph = parserGnssEph(p,navpath)
% Parse ephemeris data from RINEX file to .mat data file. 
% Supported by RINEX version 3.03
%
%%%%%-----Reference
% ftp://igs.org/pub/data/format/rinex303.pdf
% http://acc.igs.org/misc/rinex304.pdf
%%%%%-----Input
% rinex navigation file path
%
%%%%%-----Output
% Class of constellation ephemeris
%
% Author: Azurehappen
[~,~,ext] = fileparts(navpath);
if ext == ".rnx" || ext == ".obs" || ext == ".nav"
    flag = 1;
    fprintf ('Loading ephemeris...\n \n');
else
    fprintf ('File format not supported. Please input a RINEX file\n');
end
navfile = fopen(navpath);
%----------------------------------------------------%
% Initialize variables
eph.ionoParameters=[];
eph.LeapSeconds=[];
gps.prn_avb = []; gal.prn_avb = []; bds.prn_avb = []; glo.prn_avb = [];
gcount = [];ecount = [];ccount = [];rcount = [];
g = 1; e = 1; c = 1; r = 1; %Eph count
MAXGPSPRN = p.gps.max_prn; MAXGLOPRN = p.glo.max_prn;
MAXGALPRN = p.gal.max_prn; MAXBDSPRN = p.bds.max_prn;
% Create cell for t_oc, restore different size.
% The time is saved in POSIX timestamp
gps.t_oc = cell(MAXGPSPRN,1);glo.t_oc = cell(MAXGLOPRN,1);
gal.t_oc = cell(MAXGALPRN,1);bds.t_oc = cell(MAXBDSPRN,1);
% read header
fprintf ('Reading header...\n');
ionoParameters  = []; 
while (true)    
    line = fgetl(navfile);
    lineSplit = strsplit(line);
    if contains(line,'IONOSPHERIC CORR')
        if strcmp(lineSplit(1), 'GPSA')
            ionoParameters.ionoAlpha = str2double(lineSplit(2:5));
        elseif strcmp(lineSplit(1), 'GPSB')
            ionoParameters.ionoBeta = str2double(lineSplit(2:5));
        elseif strcmp(lineSplit(1), 'gal')
            ionoParameters.ionoGAL = str2double(lineSplit(2:5));
        end
            eph.ionoParameters =  ionoParameters;
    elseif contains (line,'LEAP SECONDS')
        eph.LeapSeconds = str2double(lineSplit(2));
    elseif contains(line,'END OF HEADER')
        break;
    end
end
eph.ionoParameters =  ionoParameters;
fprintf ('Finished reading the header\n \n');
%----------------------------------------------------%
% read body
fprintf ('Parsing navigation message');
while ~feof(navfile)
    line = fgetl(navfile);
    %linesp = [line(1),' ',line(2:end)]; % For example, avoid 'G15' and 'G 5' has different kind of split.
    %lineSplit = strsplit(linesp);
    sys_type = line(1);
    switch sys_type
        case{'G'} % Parse gps ephemeris data
            prn = str2double(line(2:3));
            gcount(1,g)=prn;
            gps.prn_avb(prn,1)=1; % prn_avb=1 means this satellite available in this dataset, to aviod matrix index exceed.
            indx = sum(gcount==prn);
            tget = str2double(strsplit(line(5:23)));
            %[~,~,gps.t_oc{prn}(indx)] = date2gpst(tget);
            % Time of broadcast (seconds)
            gps.t_oc{prn}(indx) = posixtime(datetime(tget));
            data = sscanf(line(24:end),'%f');
            gps.a_f0(prn,indx) = data(1); % SV clock bias (seconds) 
            gps.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec) 
            gps.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2) 
            
            data = sscanf(fgetl(navfile),'%f');  
            gps.IODE(prn,indx) = data(1); % Issue of Data, Ephemeris 
            gps.C_rs(prn,indx) = data(2); % Crs (meters) 
            gps.Delta_n(prn,indx) = data(3); % Delta n (radians/sec) 
            gps.M_0(prn,indx) = data(4); % M0 (radians) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            gps.C_uc(prn,indx) =  data(1); % Cuc (radians)
            gps.e(prn,indx) = data(2); % e Eccentricity 
            gps.C_us(prn,indx) = data(3); % Cus (radians) 
            gps.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m)) 
            
            data = sscanf(fgetl(navfile),'%f'); 
            gps.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of gps week)
            gps.C_ic(prn,indx) = data(2); % Cic (radians)
            gps.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
            gps.C_is(prn,indx) = data(4); % Cis (radians)
                
            data = sscanf(fgetl(navfile),'%f');   
            gps.i_0(prn,indx) =  data(1); % i0 (radians)
            gps.C_rc(prn,indx) = data(2); % Crc (meters) 
            gps.Omega(prn,indx) = data(3); % omega (radians)
            gps.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec) 
                
            data = sscanf(fgetl(navfile),'%f');     
            gps.IDOT(prn,indx) = data(1); % IDOT (radians/sec) 
            gps.CodesOnL2(prn,indx) = data(2); % Codes on L2 channel 
            gps.week_num(prn,indx) = data(3); % gps Week # 
            gps.L2Pflag(prn,indx) = data(4); % L2 P data flag
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            gps.SV_acc(prn,indx) = data(1); % SV accuracy (meters) See gps ICD 200H Section 20.3.3.3.1.3 
            gps.SV_health(prn,indx) = data(2); % SV health
            gps.TGD(prn,indx) = data(3); % TGD (seconds) 
            gps.IODC(prn,indx) = data(4); % IODC Issue of Data, Clock 
                
            data = sscanf(fgetl(navfile),'%f'); 
            gps.trans_time(prn,indx) = data(1); % Transmission time of message 
            if length(data) > 1
                gps.fit_interval(prn,indx) = data(2); % Fit Interval in hours
            end
            g = g +1;
            
        case{'E'} % Parse gal ephemeris data
            prn = str2double(line(2:3));
            gal.prn_avb(prn,1)=1;
            ecount(1,e)=prn;
            indx = sum(ecount==prn);
            tget = str2double(strsplit(line(5:23)));
            %[~,~,gal.t_oc{prn}(indx)] = date2gpst(tget);
            % Time of broadcast (seconds)
            gal.t_oc{prn}(indx) = posixtime(datetime(tget));
            data = sscanf(line(24:end),'%f');
            gal.a_f0(prn,indx) = data(1); % SV clock bias (seconds) 
            gal.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec) 
            gal.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2) 
            
            data = sscanf(fgetl(navfile),'%f');   
            gal.IODE(prn,indx) = data(1); % Issue of Data, Ephemeris 
            gal.C_rs(prn,indx) = data(2); % Crs (meters) 
            gal.Delta_n(prn,indx) = data(3); % Delta n (radians/sec) 
            gal.M_0(prn,indx) = data(4); % M0 (radians) 
                
            data = sscanf(fgetl(navfile),'%f'); 	  
            gal.C_uc(prn,indx) = data(1); % Cuc (radians)
            gal.e(prn,indx) = data(2); % e Eccentricity 
            gal.C_us(prn,indx) = data(3); % Cus (radians) 
            gal.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m)) 
            
            data = sscanf(fgetl(navfile),'%f'); 
            gal.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of gal week)
            gal.C_ic(prn,indx) = data(2); % Cic (radians)
            gal.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
            gal.C_is(prn,indx) = data(4); % Cis (radians)
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            gal.i_0(prn,indx) =  data(1); % i0 (radians)
            gal.C_rc(prn,indx) = data(2); % Crc (meters) 
            gal.Omega(prn,indx) = data(3); % omega (radians)
            gal.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec) 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            gal.IDOT(prn,indx) = data(1); % IDOT (radians/sec) 
            gal.Data_source(prn,indx) = data(2); % Data sources
            gal.week_num(prn,indx) = data(3); % gal Week # 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            gal.SV_acc(prn,indx) = data(1); % SISA Signal in space accuracy (meters) 
            gal.SV_health(prn,indx) = data(2); % SV health
            gal.BGD_E5a(prn,indx) = data(3); % BGD E5a/E1 (seconds) 
            gal.BGD_E5b(prn,indx) = data(4); % BGD E5b/E1 (seconds)
                
            data = sscanf(fgetl(navfile),'%f'); 
            gal.trans_time(prn,indx) = data(1); % Transmission time of message 
            e = e +1;
            
        case{'R'}
            prn = str2double(line(2:3));
            glo.prn_avb(prn,1)=1;
            rcount(1,r)=prn;
            indx = sum(rcount==prn);
            tget = str2double(strsplit(line(5:23)));
            %tget = datetime(tget); % Time of broadcast, glo time
            %tget = [tget.Year,tget.Month,tget.Day,tget.Hour,tget.Minute,tget.Second];
            %[~,~,glo.t_oc{prn}(indx)] = date2gpst(tget);
            % Represent glo time by gps time
            glo.t_oc{prn}(indx) = posixtime(datetime(tget));
            data = sscanf(line(24:end),'%f');
            glo.nTauN(prn,indx) = data(1); % SV clock bias (sec) (-TauN) 
            glo.pGammaN(prn,indx) = data(2); % SV relative frequency bias (+GammaN) 
            glo.t_of(prn,indx) = data(3); % Message frame time (tk+nd*86400) in seconds of the UTC week 
                
            data = sscanf(fgetl(navfile),'%f'); 
            glo.X(prn,indx) = data(1)*1000; % Satellite position X (m) 
            glo.Xdot(prn,indx) = data(2)*1000; % velocity X dot (m/sec) 
            glo.Xacc(prn,indx) = data(3)*1000; % X acceleration (m/sec2) 
            glo.SV_health(prn,indx) = data(4); % health (0=OK) (Bn) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            glo.Y(prn,indx) = data(1)*1000; % Satellite position Y (m) 
            glo.Ydot(prn,indx) = data(2)*1000; % velocity Y dot (m/sec)
            glo.Yacc(prn,indx) = data(3)*1000; % Y acceleration (m/sec2) 
            glo.freq(prn,indx) = data(4); % frequency number(-7...+13) (-7...+6 ICD 5.1) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            glo.Z(prn,indx) = data(1)*1000; % Satellite position Z (m) 
            glo.Zdot(prn,indx) = data(2)*1000; % velocity Z dot (m/sec) 
            glo.Zacc(prn,indx) = data(3)*1000; % Z acceleration (m/sec2) 
            glo.age(prn,indx) = data(4); % Age of operation (days)
            r = r + 1;
        case{'C'}
            prn = str2double(line(2:3));
            bds.prn_avb(prn,1)=1;
            ccount(1,c)=prn;
            indx = sum(ccount==prn);
            tget = str2double(strsplit(line(5:23)));
            tget = datetime(tget); % Time of broadcast, bds time
            tget = [tget.Year,tget.Month,tget.Day,tget.Hour,tget.Minute,tget.Second];
            %[~,~,bds.t_oc{prn}(indx)] = date2gpst(tget);
            % Represent bds time by gps time
            bds.t_oc{prn}(indx) = posixtime(datetime(tget));
            data = sscanf(line(24:end),'%f');
            bds.a_f0(prn,indx) = data(1); % SV clock bias (seconds) 
            bds.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec) 
            bds.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2) 
            
            data = sscanf(fgetl(navfile),'%f');   
            bds.AODE(prn,indx) = data(1); % Age of Data, Ephemeris  
            bds.C_rs(prn,indx) = data(2); % Crs (meters) 
            bds.Delta_n(prn,indx) = data(3); % Delta n (radians/sec) 
            bds.M_0(prn,indx) = data(4); % M0 (radians) 
                
            data = sscanf(fgetl(navfile),'%f'); 	  
            bds.C_uc(prn,indx) = data(1); % Cuc (radians)
            bds.e(prn,indx) = data(2); % e Eccentricity 
            bds.C_us(prn,indx) = data(3); % Cus (radians) 
            bds.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m)) 
            
            data = sscanf(fgetl(navfile),'%f'); 
            bds.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of bds week)
            bds.IODE(prn,indx) = mod(floor(data(1)/720),240);
            bds.C_ic(prn,indx) = data(2); % Cic (radians)
            bds.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
            bds.C_is(prn,indx) = data(4); % Cis (radians)
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            bds.i_0(prn,indx) =  data(1); % i0 (radians)
            bds.C_rc(prn,indx) = data(2); % Crc (meters) 
            bds.Omega(prn,indx) = data(3); % omega (radians)
            bds.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec) 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            bds.IDOT(prn,indx) = data(1); % IDOT (radians/sec) 
            bds.week_num(prn,indx) = data(3); % BDT Week # 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            bds.SV_acc(prn,indx) = data(1); % SV accuracy (meters) See bds ICD 200H Section 20.3.3.3.1.3 
            bds.SV_health(prn,indx) = data(2); % SV health
            bds.TGD1(prn,indx) = data(3); % TGD1   B1/B3          (seconds) 
            bds.TGD2(prn,indx) = data(4); % TGD2   B2/B3          (seconds) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            bds.trans_time(prn,indx) = data(1); % Transmission time of message 
            bds.ADOC(prn,indx) = data(2); % Age of Data Clock 
            c = c +1;       
    end
end
fclose(navfile);
if length(gps.prn_avb)<MAXGPSPRN
    gps.prn_avb = [gps.prn_avb; zeros(MAXGPSPRN-length(gps.prn_avb),1)];
end
if length(gal.prn_avb)<MAXGALPRN
    gal.prn_avb = [gal.prn_avb; zeros(MAXGALPRN-length(gal.prn_avb),1)];
end
if length(glo.prn_avb)<MAXGLOPRN
    glo.prn_avb = [glo.prn_avb; zeros(MAXGLOPRN-length(glo.prn_avb),1)];
end
if length(bds.prn_avb)<MAXBDSPRN
    bds.prn_avb = [bds.prn_avb; zeros(MAXBDSPRN-length(bds.prn_avb),1)];
end
%----------------------------------------------------%
eph.gps = gps;
eph.gal = gal;
eph.glo = glo;
eph.bds = bds;
%%Save to MAT file
% save ([fileName,'_nav.mat'], 'eph');

%----------------------------------------------------%
fprintf ('\n \nEphemeris loaded:\n \n');
fprintf ('Number of gps messages parsed: %d\n', g-1);
fprintf ('Number of GLONASS messages parsed: %d\n',r-1);
fprintf ('Number of GALileo messages parsed: %d\n',e-1);
fprintf ('Number of bds messages parsed: %d\n',c-1);

end