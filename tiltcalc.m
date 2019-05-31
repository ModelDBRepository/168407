target=3;
tiltcheck(aeon)=0;
insum=0;
for z =1:4
    for y=1:scale{z}(1)
        for x=1:scale{z}(2)
            for n =1:maxoutputs
                if connex{z,y,x,n}(1)==3
                    insum=insum+connex{z,y,x,n}(5);
                end
            end
        end
    end
end
insum=insum/(y*fanini);
for z =1:4
    for y=1:scale{z}(1)
        for x=1:scale{z}(2)
            for n =1:maxoutputs
                if connex{z,y,x,n}(1)==3 && x~=3
                    tiltcheck(aeon)=tiltcheck(aeon)+(connex{z,y,x,n}(5)-insum)*(x-3);
                end
            end
        end
    end
end

tiltcheck(aeon)=tiltcheck(aeon)/(4*y)
