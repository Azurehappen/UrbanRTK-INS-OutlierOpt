function [iono_delay, mapping_m]=ustec_iono_delay_computation(p,ustec_i,elev,az,user_t,freq)
%%%% input:
%%%% p: parameters
%%%% ustec_i: USTEC map
%%%% elev: elevation angle in radians
%%%% az: azimuth angle in radians
%%%% user_t: user gps date time

%%%% output:
%%%% iono_delay: ionospheric delay for L1 measurement computed from USTEC map

%%%% references:
%%%[1]D.W.Dodd(2007): UTILITY OF IONOSPHERE AND TROPOSPHERE MODELS FOR EXTENDING THE RANGE
%%%OF HIGH-ACCURACY GPS
%%%[2]M.B.Bakry,T.C.Wisser(1993):The FAA Wide Area Differential GPS (WADGPS) Static Ionospheric Experiment

%%% initialize
iono_delay=0;
mp_fact = 1;
%%% find suitable iono map
%%%% for debugging
iono_map_idx=[];
% load('C:\Users\ladmin\Documents\DGPS_Sirius\Mfiles_and_Data\Network_DGPS\data\035_FTP\USTEC.mat')
% ustec_i=USTEC;
iono_l=length(ustec_i);
%%%%----------------%%%%
% for ii=1:iono_l
%     %%% iono map time
%     rr=ustec_i(ii).file;
%     year_i=str2double(rr(1:4));
%     month_i=str2double(rr(5:6));
%     day_i=str2double(rr(7:8));
%     hour_i=str2double(rr(9:10));
%     min_i=str2double(rr(11:12));
%     time_i=min_i;
%     %%% time of the user
%     time_user=user_t(3);
%     
%     if (time_user-time_i)>15 || (time_user-time_i)<0
%         continue;
%     else
%         if ((day_i==user_t(1)) && (hour_i==user_t(2)))
%             iono_map_idx=ii;
%             break;
%             
%         else
%             continue;
%         end
%     end
%     
% end
for ii=1:iono_l
    %%% iono map time
    Week = ustec_i(ii).Gweek;
    DOW = ustec_i(ii).Gdow;
    SOW = ustec_i(ii).Gsow;    
    if Week == user_t.week
        if (user_t.sow - SOW < p.tec_tmax*60) && (user_t.sow - SOW >= p.tec_tmin)
            iono_map_idx = ii;
            break;
        end
    elseif user_t.week - Week == 1
        t = limit_tgps(user_t.sow - SOW);
        if t<15*60
            iono_map_idx = ii;
            break;
        end
    end
    
end
if ~isempty(iono_map_idx)
    log.tdiff = (user_t.sow - SOW)/60;
%%%% multiplying factor(vertical to slant conversion)
m=(p.Re*cos(elev)/(p.Re+p.h_iono));

%%%%%% pierce point calculation
[lat_p,lon_p]=pierce_point_calc(deg2rad(p.lat_deg),deg2rad(p.lon_deg),az,elev,p);



%%%% ionspheric delay estimation method based on interpolation
%%%% ref[2]: Appendix A
%%% take current iono map
current_iono_map=[ustec_i(iono_map_idx).TEC_table];



%%%% find columns
rr=find(current_iono_map(1,2:end)/10<lon_p);
rr1=find(current_iono_map(1,2:end)/10>lon_p);
col_idx=[rr(end)+1 rr1(1)+1];

%%% find rows
rr=find(current_iono_map(2:end,1)/10<lat_p);
rr1=find(current_iono_map(2:end,1)/10>lat_p);
row_idx=[rr(end)+1 rr1(1)+1];

%%% neighbouring points
nn_points=[current_iono_map(row_idx(1), 1) current_iono_map(1,col_idx(1));
    current_iono_map(row_idx(2), 1) current_iono_map(1,col_idx(1));
    current_iono_map(row_idx(1), 1) current_iono_map(1,col_idx(2));
    current_iono_map(row_idx(2), 1) current_iono_map(1,col_idx(2));];

nn_points=nn_points/10;

%%%% compute distance of the grid points from that pierce point
p_distance(1,1)=sqrt((lat_p-nn_points(1,1))^2+(lon_p-nn_points(1,2))^2);
p_distance(2,1)=sqrt((lat_p-nn_points(2,1))^2+(lon_p-nn_points(2,2))^2);
p_distance(3,1)=sqrt((lat_p-nn_points(3,1))^2+(lon_p-nn_points(3,2))^2);
p_distance(4,1)=sqrt((lat_p-nn_points(4,1))^2+(lon_p-nn_points(4,2))^2);


%%%% eqn (A.1)
%%%%% compute normalized inverse distance weight vector
w_dist(1,1)=1/p_distance(1,1);
w_dist(2,1)=1/p_distance(2,1);
w_dist(3,1)=1/p_distance(3,1);
w_dist(4,1)=1/p_distance(4,1);

%%%% eqn (A.1)
w_tot=sum(w_dist);

%%%% eqn (A.1)
w_dist=w_dist./w_tot;

%%%% get the respective iono vtec
vtec_npoints=[current_iono_map(row_idx(1),col_idx(1));
    current_iono_map(row_idx(2),col_idx(1));
    current_iono_map(row_idx(1),col_idx(2));
    current_iono_map(row_idx(2),col_idx(2))];

%%%% eqn (A.1)
%%% get iono vtec for that point
vtec_pp=w_dist'*vtec_npoints;

%%%% multiplying factor(vertical to slant conversion)
%%%% ref[1] (C.12.2, page-245) and ref[2] (A.3a)
m=(p.Re*cos(elev)/(p.Re+p.h_iono));
mp_fact=(1-(m)^2)^(-1/2);

%%%% get slant tec
%    stec_pp=vtec_pp/(mp_fact*10);
stec_pp=(vtec_pp/10)*mp_fact;
log.mp_fact = mp_fact;
log.vtec = vtec_pp/10;
%%% get slant iono delay
%    iono_delay=((40.3/(p.gps.L1*p.gps.L2))*stec_pp)*1e16;
%iono_delay=((40.3/(freq)^2 + 5.6*10e7/(freq)^3)*stec_pp)*1e16;
iono_delay=(40.3e16/(freq)^2)*stec_pp;
mapping_m = mp_fact * 40.3e16 / (freq * freq);
end

