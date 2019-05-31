raster=zeros(scale{2}(1),ttime);
tic
for xxxx=1:ttime
    if xxxx/10000==round(xxxx/10000)
        toc
        xxxx
    end
       for tto=1:10
        if version==2
            raster(tto,xxxx)=mem{4,xxxx}(tto,1);
        else
            raster(tto,xxxx)=mem{4,xxxx}(tto,1);
        end
       end

end
figure
space=image((raster+1)*200);title('raster')