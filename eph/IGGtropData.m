function coeff = IGGtropData(file)

nHPara = fix(6);   nPara = fix(5);
dlon = 2.5;   dlat = 2.5;
nLon = fix(360/dlon);   nLat = fix(180/dlat)+1;

fid = fopen(file);
coeff = fread(fid, 'float');
fclose(fid);
coeff = reshape(coeff, [nLon nLat nHPara nPara]);
end 