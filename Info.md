# MultiGNSS_Matlab
A MATLAB tutorial for Multi-GNSS in Standard, DGNSS, PPP

Current targeted signal:
GPS: L1
GLO:
GAL:
BDS: B1I (Before customized the satellite position calculation, please read BDS ICD, different signal channel may have different computation.)

Check BDS satellite type: http://mgex.igs.org/IGS_MGEX_Status_BDS.php

In parser_eph, the unit in GLO ephemeris related to 'km' has been convert to meter by multipling 1000.
Time in obs data need to be restricted by GPS time.

Getting USTEC data: https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/
Realtime USTEC： https://services.swpc.noaa.gov/text/us-tec-total-electron-content.txt


'.obs' and '.nav' files are in \data folder.
Getting the USTEC data, please extract the data folder to \data.
In main file, setting the name of this data folder.

Priors used in Robust Nonlinear estimation in userpos_LTS.m are
Std of prior Position = [1.414; 1.414; 1.414; 1; 0.25; 0.25]
Std of each GNSS measurement = 0.72 meters

Add the link of data source.  https://cddis.nasa.gov/archive/gnss/data/daily/2020/113/
For example: 
.20o/ is obs data for MultiGNSS, .20p is MultiGNSS eph data, 20n/ is GPS nav data.
.20d/ is COMPACT RINEX FORMAT (MultiGNSS).  20f/ is Beidou nav data.
.20g/ is GLONASS nav data. 20l/ is GAL nav data, 
.20p/ is MultiGNSS nav data. .20q/ is QZSS nav data,
.20d/ is RINEX 3 obs data, has crx compressed format.  see https://terras.gsi.go.jp/ja/crx2rnx.html
To check the the name corresponding to which data center, see http://www.igs.org/network


Debug notes: when implement PPP, some time duration has no IGS data correspond since IDOE not mathch
Checked tidx = tidx(abs(dtr)==min(abs(dtr))); in satpost_corrpsedR_singlefreq.m, it not exclude the postive value in dtr. Now fixed.

IONEX DATA from http://cddis.nasa.gov/gnss/products/ionex/
Some data have DCB information, like ehrg0620

DCB： ftp://ftp.aiub.unibe.ch/CODE/2020/
MultiGNSS DCB：ftp://igs.ign.fr/pub/igs/products/mgex/dcb/2020/

Broadcast Corrections in RTCM Version 3 Format： http://www.gnsser.com/Information/ViewDetails/535

Multi-GNSS code bias: ftp://ftp.gipp.org.cn/product/dcb/daybias/2021/

IGS data stream monitor:
https://bkgmonitor.gnssonline.eu/cgi-bin/bkgmonitor.cgi??mod=Monitoring&stream=&caster=products.igs-ip.net&moniMod=Verbose&moniTyp=All&lang=en&bDate=2020-07-11&bHour=00&bMin=00&bSec=00&eDate=2020-07-12&eHour=00&eMin=00&eSec=00

Useful PPP data source: https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/gnss_differential_code_bias_product.html