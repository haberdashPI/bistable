function actual_locs=loc_find(loc,n)
actual_locs=[];
if(n==1)
    actual_locs=loc;
else
    for i=1:n
        locs=loc+((i-1)*350);
        actual_locs=[actual_locs;locs'];
        clear locs
    end
end