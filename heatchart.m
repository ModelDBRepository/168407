heater=-1*ones(scale{2}(1),round(ttime/cycle));
heatz=zeros(scale{2}(1),round(ttime/cycle));
for xzzx=1000:ttime
     for y=1:10
         if mem{outfoc,xzzx}(y,1)>heater(y,round(xzzx/cycle))
            heater(y,round(xzzx/cycle))=mem{outfoc,xzzx}(y,1);
            if heater(y,round(xzzx/cycle))>0
                heatz(y,round(xzzx/cycle))=1;
            end
         end
     end
end
%responcestats [bandwidth,respcenter,meanstr,tilt,respsum]
if version==2
for y=1:10
    respcenter=0;
    respsum=0;
    rnum=0;
    for xzx=1:round(ttime/cycle)
        if heatz(y,xzx)==1
        
        if actdelay(xzx)>0
            updn=1;
            respcenter=respcenter+1;
        end
        if actdelay(xzx)<0
            updn=-1;
            respcenter=respcenter+1;
        end
        if actdelay(xzx)==0
            updn=0;
        end
        respsum=respsum+updn;
        rnum=rnum+1;
        end
    end
    respcenter=(respcenter);
    responcestats{aeon}(y,1)=sum(heatz(y,:));%#ok<SAGROW>
    responcestats{aeon}(y,2)=respcenter;%#ok<SAGROW>
    responcestats{aeon}(y,3)=.009+y*narrowing;%#ok<SAGROW>
    responcestats{aeon}(y,4)=tilt; %#ok<SAGROW>
    responcestats{aeon}(y,5)=respsum; %#ok<SAGROW>
end
end
if version==2
  
responcestats{aeon}
['aeon: ',aeon]
end






%figure
%subplot(1,1,1)
%heatpic=image((heater+1)*200);
%hgsave(['heat','aeon',num2str(aeon)])
