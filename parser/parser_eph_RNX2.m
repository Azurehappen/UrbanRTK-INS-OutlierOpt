function eph = parser_eph_RNX2(p,navpath)
% Parse ephemeris data from .20N file to matlab data file. 
% Supported by RINEX version 2.11
%
%%%%%-----Reference
% https://kb.igs.org/hc/en-us/articles/115003980188-RINEX-2-11
%
%%%%%-----Input
% .nav file path
%
%%%%%-----Output
% Class of constellation ephemeris
%
% Author: Wang Hu
% [fileName, filePath] = uigetfile({'*.nav', 'Navigation Files (*.nav)'}, 'Select a Navigation File');
% navpath = strcat (filePath, fileName);
[~,~,ext] = fileparts(navpath);

if strcmp(ext,'.20n')||strcmp(ext,'.20N')
fprintf ('Loading ephemeris...\n \n');
EndOfHeader = 0;
navfile = fopen(navpath);
%----------------------------------------------------%
% Initialize variables
eph.DOcreate=[];
eph.ionosphericParameters=[];
eph.LeapSeconds=18;
gps.svid_avb = [];
gcount = [];
g = 1; %Eph count
% read header
fprintf ('Reading header...\n');
while (~EndOfHeader)
    ionosphericParameters  = [];
    
    line = fgetl(navfile);
    lineSplit = strsplit(line);
    
    if contains(line,'RINEX VERSION')
        Version = lineSplit(2);
        if ~strcmp(Version{1}(1),'2')
            error('Not the correct version, should be 2')
        end  
%     elseif contains(line,'DATE')
%         date = lineSplit(4);
%         year = str2double(date{1,1}(1:4));
%         month = str2double(date{1,1}(5:6));
%         day = str2double(date{1,1}(7:8));
%         eph.DOcreate = [year,month,day]; % Date of file creation 
%     elseif contains(line,'IONOSPHERIC CORR')
%         if strcmp(lineSplit(1), 'GPSA')
%             ionosphericParameters.ionoAlpha = str2double(lineSplit(2:5));
%         elseif strcmp(lineSplit(1), 'GPSB')
%             ionosphericParameters.ionoBeta = str2double(lineSplit(2:5));
%         elseif strcmp(lineSplit(1), 'GAL')
%             ionosphericParameters.ionoGAL = str2double(lineSplit(2:5));
%         end
%             eph.ionosphericParameters =  ionosphericParameters;
%     elseif contains (line,'LEAP SECONDS')
%         eph.LeapSeconds = str2double(lineSplit(2));
    elseif contains(line,'END OF HEADER')
       EndOfHeader = 1;
    end
end
fprintf ('Finished reading the header\n \n');
%----------------------------------------------------%
% read body
fprintf ('Parsing navigation message');
while ~feof(navfile)
    line = fgetl(navfile);
%     switch sat_id
%         case{'G'} % Parse GPS ephemeris data
            svid = str2double(line(1:2));
            gcount(1,g)=svid;
            gps.svid_avb(svid,1)=1; %gps.svid=1 means this satellite available in this dataset, to aviod matrix index exceed.
            indx = sum(gcount==svid);
            date = strsplit(line(4:22));
            date{1} = ['20' date{1}];
            [~,~,gps.t_oc(svid,indx)] = date2gnsst(str2double(date)); % Time of broadcast (seconds)
            gps.a_f0(svid,indx) = str2double(line(23:41)); % SV clock bias (seconds) 
            gps.a_f1(svid,indx) = str2double(line(42:60)); % SV clock drift (sec/sec) 
            gps.a_f2(svid,indx) = str2double(line(61:79)); % SV clock drift rate (sec/sec2) 
            
            lineSplit = fgetl(navfile);  
            gps.IODE(svid,indx) = str2double(lineSplit(4:22)); % Issue of Data, Ephemeris 
            gps.C_rs(svid,indx) = str2double(lineSplit(23:41)); % Crs (meters) 
            gps.Delta_n(svid,indx) = str2double(lineSplit(42:60)); % Delta n (radians/sec) 
            gps.M_0(svid,indx) = str2double(lineSplit(61:79)); % M0 (radians) 
                
            lineSplit = fgetl(navfile);	  
            gps.C_uc(svid,indx) = str2double(lineSplit(4:22)); % Cuc (radians)
            gps.e(svid,indx) = str2double(lineSplit(23:41)); % e Eccentricity 
            gps.C_us(svid,indx) = str2double(lineSplit(42:60)); % Cus (radians) 
            gps.sqrtA(svid,indx) = str2double(lineSplit(61:79)); % sqrt(A) (sqrt(m)) 
            
            lineSplit = fgetl(navfile);
            gps.t_oe(svid,indx) = str2double(lineSplit(4:22)); % Toe Time of Ephemeris (sec of GPS week)
            gps.C_ic(svid,indx) = str2double(lineSplit(23:41)); % Cic (radians)
            gps.Omega_0(svid,indx) = str2double(lineSplit(42:60)); % OMEGA0 (radians)
            gps.C_is(svid,indx) = str2double(lineSplit(61:79)); % Cis (radians)
                
            lineSplit = fgetl(navfile);	    
            gps.i_0(svid,indx) =  str2double(lineSplit(4:22)); % i0 (radians)
            gps.C_rc(svid,indx) = str2double(lineSplit(23:41)); % Crc (meters) 
            gps.Omega(svid,indx) = str2double(lineSplit(42:60)); % omega (radians)
            gps.OmegaDot(svid,indx) = str2double(lineSplit(61:79)); % OMEGA DOT (radians/sec) 
                
            lineSplit = fgetl(navfile);	    
            gps.IDOT(svid,indx) = str2double(lineSplit(4:22)); % IDOT (radians/sec) 
            gps.CodesOnL2(svid,indx) = str2double(lineSplit(23:41)); % Codes on L2 channel 
            gps.week_num(svid,indx) = str2double(lineSplit(42:60)); % GPS Week # 
            gps.L2Pflag(svid,indx) = str2double(lineSplit(61:79)); % L2 P data flag
                
            lineSplit = fgetl(navfile);	    
            gps.SV_acc(svid,indx) = str2double(lineSplit(4:22)); % SV accuracy (meters) See GPS ICD 200H Section 20.3.3.3.1.3 
            gps.SV_health(svid,indx) = str2double(lineSplit(23:41)); % SV health
            gps.TGD(svid,indx) = str2double(lineSplit(42:60)); % TGD (seconds) 
            gps.IODC(svid,indx) = str2double(lineSplit(61:79)); % IODC Issue of Data, Clock 
                
            lineSplit = fgetl(navfile);
            gps.trans_time(svid,indx) = str2double(lineSplit(4:22)); % Transmission time of message 
            gps.fit_interval(svid,indx) = str2double(lineSplit(23:41)); % Fit Interval in hours
            g = g +1;
end
fclose(navfile);
if length(gps.svid_avb)<p.gps.num_prn
    gps.svid_avb = [gps.svid_avb; zeros(33-length(gps.svid_avb),1)];
end
%----------------------------------------------------%
eph.gps = gps;
%%Save to MAT file
% save ([fileName,'_nav.mat'], 'eph');
else
fprintf ('File format not suppoerted. Please input a RINEX navigation (.nav) file\n');
end
%----------------------------------------------------%
fprintf ('\n \n Ephemeris loaded correctly\n \n');
fprintf ('Number of GPS messages parsed: %d\n', g-1);

end