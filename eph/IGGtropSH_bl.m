function ZTD = IGGtropSH_bl(siteLon, siteLat, h, doy)
% (c) State Key Laboratory  of Geodesy and Earth's Dynamics, Institute of Geodesy and Geophysics, 
% Chinese Academy of Sciences 2017

% This program calculates ZTD from IGGtropSH model, which is established by Li Wei. 
% The coefficients of IGGtropSH are derived using ERA-Interim reanalysis pressure level data during 2006.1-2009.12

% References:
% Li W, Yuan YB, Ou JK, He YJ (2018). IGGtrop_SH and IGGtrop_rH: Two improved empirical tropospheric delay models based on vertical 
% reduction functions. IEEE Transactions on Geoscience and Remote Sensing. PP. 1-13. 10.1109/TGRS.2018.2812850.
% Li W, Yuan YB, Ou JK, Li H, Li ZS (2012). A new global zenith tropospheric delay model IGGtrop for GNSS applications. Chin Sci Bull, 57(17): 2132¨C2139

%==========================input============================
% sitelon:longitude [degree]; range: 0 to 360 degree
% sitelat:latitude [degree] ; range: -90 to 90 degree
% siteh:orthometric height [km]; 
% doy: day of year
%===========================================================
%==========================output===========================
% ZTD: zenith tropospheric delay [m]
%===========================================================

% example
% siteLon = 128.9;
% siteLat = -90;
% h = 0.232;
% doy = 1:1:365;

nHPara = fix(6);   nPara = fix(5);
dlon = 2.5;   dlat = 2.5;
nLon = fix(360/dlon);   nLat = fix(180/dlat)+1;
lon = dlon*(0:1:nLon-1);   lat = 90-dlat*(0:1:nLat-1);

fid = fopen('IGGtropSHexpModel.ztd');
coeff = fread(fid, 'float');
fclose(fid);
coeff = reshape(coeff, [nLon nLat nHPara nPara]);

ix = (siteLon-lon(1))/(lon(2)-lon(1)) + 1;
iy = (siteLat-lat(1))/(lat(2)-lat(1)) + 1;

Y_grid = zeros( size(doy,2), 4);

for igrid=1:4
    a=mod( fix(ix+bitand(igrid-1, 1))-1, nLon) + 1;
    b=fix(iy)+ fix((igrid-1)/2);
    if b>nLat
        b=nLat;
    end 
    a_coeff = coeff(fix(a), fix(b), :, :);
    a_coeff = reshape(a_coeff, [nHPara nPara]);
    
    a=(a_coeff(2:nHPara-1, 1))';
    b=h.^(1:1:nHPara-2)';
    ave=a_coeff(1,1)*exp(a*b) + a_coeff(fix(nHPara),1);

    a=a_coeff(:,2:nPara)';
    b=h.^(0:1:nHPara-1)';
    A=a*b;    
    t=2*pi/365.25*doy;
    Y_grid(:, igrid) = ave(1)+A(1)*cos(t)+A(2)*sin(t)+A(3)*cos(2*t)+A(4)*sin(2*t);    
end
p = ix-fix(ix);
q = iy-fix(iy);
ZTD = (1-p)*(1-q)*Y_grid(:,1) + (1-q)*(p)*Y_grid(:,2) +  (1-p)*(q)*Y_grid(:,3)+ (p)*(q)*Y_grid(:,4);






