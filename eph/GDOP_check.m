function tab = GDOP_check(p,GDOP)
    tab = 1;
    if GDOP>p.GDOP_mark
        tab = 0;
    end
end