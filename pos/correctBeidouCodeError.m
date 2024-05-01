function dataStruct = correctBeidouCodeError(p,dataStruct)

for i = 1:length(dataStruct.svprn_mark)
    if dataStruct.svprn_mark(i) == p.bds.sys_num
        corr = beidouCodeCorr(dataStruct.prn_record(i), ...
            rad2deg(dataStruct.elev(i)));
        dataStruct.corr_range(i) = dataStruct.corr_range(i) + corr;
    end
end

end