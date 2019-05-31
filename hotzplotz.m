svbt=zeros(10,2);
figure
hold on    

for aeon=1:2:aeons
    for y=1:10
        svbt(y,1)=svbt(y,1)+responcestats{aeon}(y,3);
        svbt(y,2)=svbt(y,2)+responcestats{aeon}(y,1);
        plot(responcestats{aeon}(y,3),responcestats{aeon}(y,1))
    end

end
svbt=svbt/aeons;
plot(svbt(:,1),svbt(:,2))

tvmt=zeros(aeons/2,3);
tvm=zeros(aeons/2,3);
    figure
    hold on

for y=1:10

    for aeon=1:2:aeons
        
        tvm((aeon+1)/2,1)=responcestats{aeon}(y,4);
        tvm((aeon+1)/2,2)=responcestats{aeon}(y,2);
        tvm((aeon+1)/2,3)=responcestats{aeon}(y,5);
    end
    tvmt=tvmt+tvm;
end
tvmt=tvmt/y;
dsi=tvmt(:,3)./tvmt(:,2);
plot(dsi)
figure
hold on
%plot(-.0015:.0001:.0015,dsi)
plot(dsi)
plot(34.75*tiltcheck(1:2:end),'r')
hgsave([num2str(aeon),'dsivtilt'])
