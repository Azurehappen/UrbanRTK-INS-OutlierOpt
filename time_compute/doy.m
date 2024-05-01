function num = doy(date)
% Compute day of year
% date: [year,month,day]
leap_year = [31,29,31,30,31,30,31,31,30,31,30,31];
common = [31,28,31,30,31,30,31,31,30,31,30,31];
if (mod(date(1),4)==0 && mod(date(1),100)~=0)|| mod(date(1),400)==0
    num = sum(leap_year(1:date(2)-1))+date(3);
else
    num = sum(common(1:date(2)-1))+date(3);
end

end