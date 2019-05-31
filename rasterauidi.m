
for xxxx=1:ctime
    for yyyy=1:scale{2}(1)
        raster(yyyy,xxxx)=mem{2,xxxx}(yyyy);
    end
end
space=image((raster+.94)*200);title('raster')