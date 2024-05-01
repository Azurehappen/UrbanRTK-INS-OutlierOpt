function bia = parser_bia(p,biapath)
% Parse code bias data from OSR bia file to matlab data file.
% GPS/GLO/BDS system can also be obtained from CNES SSR data
biafile = fopen(biapath);
while (true)
    line = fgetl(biafile);
    if contains(line,'+BIAS/SOLUTION')
        line = fgetl(biafile);
        break;
    end
end
PRN_i = strfind(line,'PRN');
type_i = strfind(line,'OBS1');
val_i = strfind(line,'VALUE____');
bia.GPS.bia_C1C = NaN(p.gps.max_prn,1);
bia.GPS.bia_C2L = NaN(p.gps.max_prn,1);
bia.GAL.bia_C1C = NaN(p.gal.max_prn,1);
bia.GAL.bia_C7Q = NaN(p.gal.max_prn,1);
bia.BDS.bia_C2I = NaN(p.bds.max_prn,1);
bia.BDS.bia_C7I = NaN(p.bds.max_prn,1);
while ~feof(biafile)
    line = fgetl(biafile);
    if strcmp(line(PRN_i),'G')
        prn = str2double(line(PRN_i+1:PRN_i+2));
        type = line(type_i:type_i+3);
        if contains(type,'C1C')
            bia.GPS.bia_C1C(prn,1) = str2double(line(val_i:val_i+8));
        end
        if contains(type,'C2L')
            bia.GPS.bia_C2L(prn,1) = str2double(line(val_i:val_i+8));
        end
    end
    if strcmp(line(PRN_i),'E')
        prn = str2double(line(PRN_i+1:PRN_i+2));
        type = line(type_i:type_i+3);
        if contains(type,'C1C')
            bia.GAL.bia_C1C(prn,1) = str2double(line(val_i:val_i+8));
        end
        if contains(type,'C7Q')
            bia.GAL.bia_C7Q(prn,1) = str2double(line(val_i:val_i+8));
        end
    end
    if strcmp(line(PRN_i),'C')
        prn = str2double(line(PRN_i+1:PRN_i+2));
        type = line(type_i:type_i+3);
        if contains(type,'C2I')
            bia.BDS.bia_C2I(prn,1) = str2double(line(val_i:val_i+8));
        end
        if contains(type,'C7I') && prn <= 18
            bia.BDS.bia_C7I(prn,1) = str2double(line(val_i:val_i+8));
        end
        if contains(type,'C7Z') && prn > 18
            bia.BDS.bia_C7I(prn,1) = str2double(line(val_i:val_i+8));
        end
    end
end
fclose(biafile);