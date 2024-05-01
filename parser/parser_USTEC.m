function USTEC = parser_USTEC(USTEC_path)
% ------------------------
% Input: 
% folder of the USTEC files
% 
% Description: 
% This script is for parsing the USTEC VTEC maps extracted from zipped files
% downloaded from the archive of USTEC: https://www.ngdc.noaa.gov/stp/IONO/USTEC/products/
% User can extract the VTEC map files to a folder and use the folder as the
% input.
% 
% How to use: Base on the name of the folder and files, user may need to
%             change the address of the utility folder and input file folder.
% --Zeyi Jiang & Wang Hu
% ------------------------


USTEC_file = dir(fullfile(USTEC_path,'*_TEC.txt'));
USTEC=struct;

for i_file = 1:length(USTEC_file)
    filename_USTEC = strcat(USTEC_file(i_file).folder,'\',USTEC_file(i_file).name);
    fileID = fopen(filename_USTEC);
    flag_log = 0;
    start_line = 0;
    line_now = 0;
    while ~feof(fileID)
%         text_lines = text_lines +1;
        line_now = line_now+1;
        data_line = fgetl(fileID);
      
        if contains(data_line,':Product: ')
            split_cell = strsplit(data_line);
            time = split_cell{2};
            year = str2double(time(1:4));
            mon = str2double(time(5:6));
            day = str2double(time(7:8));
            hour = str2double(time(9:10));
            min = str2double(time(11:12));
            sec =0;
            [gps_week, gps_dow, gps_sow] = utc2gpst([year,mon,day,hour,min,sec]);
        end
        
        %%%Find the begining 
        if contains(data_line,'#                 Vertical and Slant Path Total Electron Content')
           start_line = line_now+3;
        end
        if line_now == start_line
           flag_log = 1; 
        end
        %%%logging 
        if flag_log == 1
            if strcmp(string(data_line(1:3)),'999')
                break
            end
            split_cell = strsplit(data_line);
            num=length(split_cell);
            for ii=1:num
                thisline(ii)=str2double(split_cell{1,ii});
            end
            thisline = thisline(2:end-1);
            try
                USTEC(i_file).TEC_table = [USTEC(i_file).TEC_table;thisline];
                USTEC(i_file).file = USTEC_file(i_file).name;
            catch
                USTEC(i_file).TEC_table = thisline;
                USTEC(i_file).file = USTEC_file(i_file).name;
                USTEC(i_file).Gweek = gps_week;
                USTEC(i_file).Gdow = gps_dow;
                USTEC(i_file).Gsow = gps_sow;
            end
        end
        
        
    end
    fclose(fileID);
end

end