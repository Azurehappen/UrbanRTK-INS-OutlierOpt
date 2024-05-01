function load_USTEC(st,et,load_path)
% ------------------------
% Input: 
% Start time and end time
% 
% Description: 
% This script is for parsing the USTEC VTEC from 
% https://www.ngdc.noaa.gov/
%
% --Zeyi Jiang & Wang Hu
% ------------------------
fprintf ('Downloading USTEC data...\n');
weboptions.Timeout = 30; % Timeout 30s
st.Minute = floor(st.Minute/15)*15;
st.Second = 0;
et.Minute = ceil(et.Minute/15)*15;
et.Second = 0;
if st.Day ~= et.Day || days(et-st)>1
    mid = st; mid = mid + days(1);
    if mid.Hour ~= 0 || mid.Minute ~= 0 || mid.Second ~= 0
        mid.Hour = 0; mid.Minute = 0; mid.Second = 0;
    end
    t_sr = mid:days(1):et;
    if mid ~= st
        t_sr = [st,t_sr];
    end
    if t_sr(end)~=et
        t_sr = [t_sr,et];
    end
    for i = 1:length(t_sr)-1
        t = t_sr(i);
        Year = string(t,'yyyy'); Mon = string(t,'MM');
        Day = string(t,'dd');
        url = "https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/" + Year + "/"...
                + Mon + "/" + Day;
        % Check if url exist
        try
            S = webread(url);
        catch
            S = -1;
        end
        if S == -1
            % Download compressed data file
            url = "https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/" + Year + "/"...
                + Mon + "/" + string(t,'yyyyMMdd') + "_ustec.tar.gz";
            file = load_path + string(t,'yyyyMMdd') + "_ustec.tar.gz";
            gzfile = websave(file,url);
            tarfile = gunzip(gzfile);
            delete(gzfile);
            untar(tarfile);
            delete(tarfile);
        else
            % Download data file one by one
            for tt = t:minutes(15):t_sr(i+1)-minutes(15)
                Year = string(tt,'yyyy'); Mon = string(tt,'MM');
                Day = string(tt,'dd');
                url = "https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/" + Year + "/"...
                        + Mon + "/" + Day + "/" + string(t,'yyyyMMddHHmm') + "_TEC.txt";
                file = load_path + string(t,'yyyyMMddHHmm') + "_TEC.txt";
                websave(file,url);
            end
        end
    end
else
    Year = string(st,'yyyy'); Mon = string(st,'MM');
    Day = string(st,'dd');
    url = "https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/" + Year + "/"...
        + Mon + "/" + Day;
    % Check if url exist
    try
        S = webread(url);
    catch
        S = -1;
    end
    if S == -1
        % Download compressed data file
        url = "https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/" + Year + "/"...
            + Mon + "/" + string(st,'yyyyMMdd') + "_ustec.tar.gz";
        file = load_path + string(st,'yyyyMMdd') + "_ustec.tar.gz";
        gzfile = websave(file,url);
        tarfile = gunzip(gzfile);
        delete(gzfile);
        untar(tarfile);
        delete(tarfile);
    else
        for t = st:minutes(15):et
            Year = string(t,'yyyy'); Mon = string(t,'MM');
            Day = string(t,'dd');
            url = "https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/" + Year + "/"...
                + Mon + "/" + Day + "/" + string(t,'yyyyMMddHHmm') + "_TEC.txt";
            file = load_path + string(t,'yyyyMMddHHmm') + "_TEC.txt";
            if exist(file,'file')~=2
                websave(file,url);
            end
        end
    end

end
fprintf ('Down.\n');
