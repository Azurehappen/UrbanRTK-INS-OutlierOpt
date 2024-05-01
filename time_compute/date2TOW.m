function iTOW = date2TOW(y,m,d,h,min,s)
    day=strcat(string(2000+y),'-',string(m),'-',string(d));
    [DayNum,~]=weekday(day);
    iTOW = (DayNum-1)*24*3600+h*3600+min*60+s;
end