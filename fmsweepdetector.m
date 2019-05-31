%1. initialization

clear %clears existing variables
%======= version selects which FMsweep detector is being used.
version=4; % 1==facilitation 2==sideband 3==test 4==duration tuning
shutter=1;opti=0;
tic %starts timer

graph=1; %decides if we create a graph
moviem=0; % dectides if we make a movie
if version==1
    scale{1}=[10,10];
    scale{2}=[10,10];
    scale{3}=[10,5];
    scale{4}=[10,1];
end
if version==2
    scale{1}=[1,5];%y=1
    scale{2}=[10,1];
    scale{3}=[10,1];
    scale{4}=[10,1];
end
if version==3
    scale{1}=[10,1];
    scale{2}=[10,1];
    scale{3}=[10,1];
    scale{4}=[10,1];
end
if version==4
    scale{1}=[10,5];
    scale{2}=[10,5];
    scale{3}=[10,5];
    scale{4}=[10,5];
    %scale{1}=[10,5];
    %scale{2}=[10,5];
    %scale{3}=[10,5];
    %scale{4}=[10,5];
end
stdp=0;
ostdp=stdp;
training=0;
balance_input=1;
%training_reliability=.6;
window=100000;
lag(1)=10; %timing of inputs for second even
lag(2)=10; %timing of inputs for third event
lag(3)=10;
lag(4)=10;
lag(5)=-10;
tfires=0;
events=5;
cycle=2000;
lowest=0;

preset=lowest-100;
failed=0; %error checking
%maxsyns=50; %max numbe of synapic events a cell can process simultinously
maxoutputs = 4000; %number of outgoing synapses
lmax=4; % number of layers

% set size of layers
if version==2
    for z=1:4
        syncount{z}=zeros(scale{z}(1),scale{z}(2));
    end
end


tau=40;


al = 3.65;  %map stretch
sig = 0.06; %baseline input
mu = 0.0005; %scales change of u
beta_e = 0.133;% scales input
sscale=0.001;
wscale=-0.002;
aeons=1;

if version<=2
    stackmax=2000;
    for ii=1:stackmax
        countdown(ii,:)=[0 0 0 1 0];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. connectivity (CPG)
for z=1:lmax
    for y=1:scale{z}(1)
        for x=1:scale{z}(2)
            synsum(z,y,x)=0;
            for n=1:maxoutputs
                connex{z,y,x,n}=zeros(1,5); % % set up connection matrix 25 connections per cell, 3 connectivity dimensions z,y,x
            end
        end
    end
end

if version==4
    events=9;
    tticks=0;
    for ducks=1:events
        for tics=1:25:((25*ducks))%1:25:(125+(25*ducks))
            tticks=tticks+1;
            event{tticks}=[((ducks*cycle)+tics),(ducks-8)];
        end
        lowest=0;
        
    end
    events=tticks;
    %layer1==>layer2
    for y=1:scale{2}(1)
        for x=1:1:scale{2}(1)
            if y==1
                connex{1,y,x,1}=[2,y,x,.01,.012]; %#ok<*SAGROW> %.012
            end
            if y==2
                connex{1,y,x,1}=[2,y,x,.01,.0115]; %#ok<*SAGROW> %.010
            end
            if y==3
                connex{1,y,x,1}=[2,y,x,.01,.010]; %#ok<*SAGROW> %.008
            end
            if y==4
                connex{1,y,x,1}=[2,y,x,.01,.0087]; %#ok<*SAGROW> %.007
            end
            if y==5
                connex{1,y,x,1}=[2,y,x,.01,.008]; %#ok<*SAGROW> %.0068
            end
            if y==6
                connex{1,y,x,1}=[2,y,x,.01,.00745]; %#ok<*SAGROW> %.012
            end
            if y==7
                connex{1,y,x,1}=[2,y,x,.01,.007]; %#ok<*SAGROW> %.010
            end
            if y==8
                connex{1,y,x,1}=[2,y,x,.01,.00695]; %#ok<*SAGROW> %.008
            end
            if y==9
                connex{1,y,x,1}=[2,y,x,.01,.0068]; %#ok<*SAGROW> %.007
            end
            if y==10
                connex{1,y,x,1}=[2,y,x,.01,.0050]; %#ok<*SAGROW> %.0068
            end
        end
    end
    
    
    %layer1==>layer3
    for y=1:scale{2}(1)
        for x=1:scale{1}(2)
            connex{1,y,x,2}=[3,y,x,.01,.085];
        end
    end
    
    %layer2==>layer4
    for y=1:scale{2}(1)
        for x=1:scale{1}(2)
            connex{2,y,x,1}=[4,y,x,.01,.06];%0.06
        end
    end
    
    %layer3==>layer4
    for y=1:scale{2}(1)
        for x=1:scale{3}(2)
            connex{3,y,x,1}=[4,y,x,.01,-.033];%-.033
        end
    end
end
if version==3
    
    for ducks=1:events
        event{ducks}=[ducks*cycle,(-4*ducks)+4];%(ducks-16)]; % when event happens
        if (event{ducks}(2))*scale{1}(2)<lowest
            lowest=(event{ducks}(2))*scale{1}(2);
        end
    end
    connex{1,1,1,1}=[2,1,1,1,.03];
    connex{2,1,1,1}=[3,1,1,1,-.2];
    connex{3,1,1,1}=[4,1,1,1,.1];
end
if version==2
    faninx=5;
    fanini=5;
    latsa=5;
    %layer1==>layer2
    z=1;
    for y=1:scale{2}(1)
        for x=1:scale{1}(2)
            for n=0:(faninx-1)
                if x-n+((faninx-1)/2)-((scale{1}(2)-scale{2}(2))/2)>0 && x-n+((faninx-1)/2)-((scale{1}(2)-scale{2}(2))/2)<=(scale{2}(2))
                    syncount{z}(1,x)=syncount{z}(1,x)+1;
                    connex{1,1,x,syncount{z}(1,x)}=[2,y,x-n+((faninx-1)/2)-((scale{1}(2)-scale{2}(2))/2),5,0.0075];%.01%.012 %#ok<SAGROW>
                end
            end
        end
    end
    %layer2==>layer4
    z=2;
    for y=1:scale{2}(1)
        for x=1:scale{2}(2)
            syncount{z}(y,x)=syncount{z}(y,x)+1;
            connex{2,y,x,syncount{z}(y,x)}=[4,y,1,5,.0035]; %#ok<SAGROW>
        end
    end
    %     %layer2==>layer2
    %     for y=1:scale{2}(1)
    %         for x=1:scale{2}(2)
    %             for n=1:latsa
    %                 if x+n-3>0 && x+n-3<=scale{2}(2)
    %                     connex{2,y,x,n+1}=[2,y,x+n-3,5,.002];
    %                     %[2,y,x,n+1]
    %                     %connex{2,y,x,n+1}
    %                 end
    %                 if x+n-3>0 && x+n-3<=scale{3}(2)
    %                     connex{2,y,x,n+1+latsa}=[3,y,x+n-3,5,.002];
    %                 end
    %             end
    %         end
    %     end
    % %layer1==>layer3
    z=1;
    narrowing=.0005;
    tilt=0%0.0001*(aeon-16);%0.00%1%
    for y=1:scale{3}(1)
        for x=1:scale{1}(2)
            for n=1:fanini
                if (x+n-fanini)>0 && (x+n-fanini)<(scale{3}(2)+1)
                    syncount{z}(1,x)=syncount{z}(1,x)+1;
                    connex{1,1,x,syncount{z}(1,x)}=[3,y,x+n-fanini,5,.009+y*narrowing+(n-(fanini/2))*tilt]; %#ok<SAGROW>
                end
            end
        end
    end
    
    
    
    %layer3==>layer2
    z=3;
    for y=1:scale{2}(1)
        for x=1:scale{3}(2)
            if x>0 && x<=scale{2}(2)
                syncount{z}(y,x)=syncount{z}(y,x)+1;
                connex{3,y,x,syncount{z}(y,x)}=[2,y,x,5,-0.2];
            end
        end
    end
end
if version==1
    events=20;
    for ducks=1:events
        event{ducks}=[ducks*cycle,(-1*ducks)+14];%(ducks-16)]; % when event happens
        if (event{ducks}(2))*scale{1}(2)<lowest
            lowest=(event{ducks}(2))*scale{1}(2);
        end
    end
    %layer1==>layer2
    for y=1:scale{2}(1)
        for x=1:scale{1}(2)
            if x-1>0
                
                connex{1,1,x,(2+(y-1)*2)}=[2,y,x-1,y,.00806];%.00806
            end
            if x<=scale{2}(2)
                connex{1,1,x,(1+(y-1)*2)}=[2,y,x,10,.008057];%.008057
            end
        end
    end
    
    % %layer1==>layer3
    %for x=1:scale{1}(2)
    %end
    
    %layer2==>layer4
    for y=1:scale{2}(1)
        for x=1:scale{2}(2)
            connex{2,y,x,1}=[4,y,1,50,.0025];
        end
    end
    
    %layer3==>layer2
    for x=1:scale{3}(2)
        for y=1
            if x>0 && x<=scale{2}(2)
                connex{3,y,x,n}=[2,y,x,50,0];
            end
        end
    end
end
for z=1:lmax %for all layers
    for y= 1:scale{z}(1)
        for x= 1:scale{z}(2)
            startsynsum(z,y,x)=0;
        end
    end
end
if balance_input==1
    for z=1:lmax %for all layers
        for y= 1:scale{z}(1)
            for x= 1:scale{z}(2)
                for n=1:maxoutputs
                    k=connex{z,y,x,n}(1);
                    j=connex{z,y,x,n}(2);
                    i=connex{z,y,x,n}(3);
                    s=connex{z,y,x,n}(4);
                    p=connex{z,y,x,n}(5);
                    if k~=0 && p>0
                        startsynsum(k,j,i)=startsynsum(k,j,i)+p;
                    end
                end
            end
        end
    end
end






for aeon=1:aeons
    close all
    if version==2
        clear event
        if training==0 || aeon/2~=round(aeon/2)
            events=141;%70
            for ducks=1:events
                actdelay(ducks)=(1*ducks)-71;
                event{ducks}=[ducks*cycle,actdelay(ducks)];%(ducks-16)]; % when event happens 36
                if (event{ducks}(2))*scale{1}(2)<lowest
                    lowest=(event{ducks}(2))*scale{1}(2);
                end
            end
            if training==1;
                stdp=0;
            end
        else
            stdp=ostdp;
            events=10;
            for ducks=1:events
                if ducks<=events*training_reliability
                    event{ducks}=[ducks*cycle,(-8)];%(ducks-16)]; % when event happens 36
                    if (event{ducks}(2))*scale{1}(2)<lowest
                        lowest=(event{ducks}(2))*scale{1}(2);
                        
                    end
                else
                    event{ducks}=[ducks*cycle,(8)];%(ducks-16)]; % when event happens 36
                    if (event{ducks}(2))*scale{1}(2)<lowest
                        lowest=(event{ducks}(2))*scale{1}(2);
                        
                    end
                end
            end
        end
    end
    ttime=event{events}(1)+400; % number of iterations
    view=100000;
    eye=zeros(1,ttime);
    blink=50;
    
    for zz=1:lmax % going thought all layers
        for ctime=1:ttime %for all time
            mem{zz,ctime}=zeros(scale{zz});  % pre-alocates memory matrix
        end
        cellset
        input{zz}=zeros(scale{zz});
        if stdp==1
            lastfired{zz,scale{zz}(1),scale{zz}(2)}=0;
            
        end
        
        fcount{zz}=zeros(scale{zz});
        for cevent=1:5
            fcounts{zz,cevent}=zeros(scale{zz});
        end
        if version<=2
            delaystack{zz}=zeros(stackmax,5);
        end
        if opti==0
            coef{zz}=0.2*ones(scale{zz});
        end
    end
    evec=ones((events*scale{1}(2))+1,2)*ttime+1;
    ticker=1;
    empt=0;
    for zxy=1:events
        for xzy=1:scale{1}(2)
            evec(ticker,1)=(event{zxy}(1)+event{zxy}(2)*xzy);
            evec(ticker,2)=xzy;
            ticker=ticker+1;
        end
    end
    bean=1;
    evec=sortrows(evec);
    ttime %#ok<NOPTS>
    for ctime=3:ttime %for almost all time
        if ctime/1000==(round(ctime/1000))
            ctime %#ok<NOPTS>
            toc
            %'run speed' %#ok<NOPTS>
            %ctime/toc %#ok<NOPTS>
        end
        
        %4(a) reset factors
        if tfires>0
            for z=1:lmax %for all layers
                for y= 1:scale{z}(1)
                    for x= 1:scale{z}(2)
                        synsum(z,y,x)=0;
                    end
                end
            end
            if balance_input==1
                for z=1:lmax %for all layers
                    for y= 1:scale{z}(1)
                        for x= 1:scale{z}(2)
                            for n=1:maxoutputs
                                k=connex{z,y,x,n}(1);
                                j=connex{z,y,x,n}(2);
                                i=connex{z,y,x,n}(3);
                                s=connex{z,y,x,n}(4);
                                p=connex{z,y,x,n}(5);
                                if k~=0 && p>0
                                    synsum(k,j,i)=synsum(k,j,i)+p;
                                end
                            end
                        end
                    end
                end
            end
            if balance_input==1
                for z=1:lmax %for all layers
                    for y= 1:scale{z}(1)
                        for x= 1:scale{z}(2)
                            for n=1:maxoutputs
                                k=connex{z,y,x,n}(1);
                                j=connex{z,y,x,n}(2);
                                i=connex{z,y,x,n}(3);
                                s=connex{z,y,x,n}(4);
                                p=connex{z,y,x,n}(5);
                                if k~=0 &&p>0
                                    connex{z,y,x,n}(5)=connex{z,y,x,n}(5)*(startsynsum(k,j,i)/synsum(k,j,i));
                                end
                            end
                        end
                    end
                end
            end
        end
        tfires=0;
        for z=1:lmax
            
            %this block takes the current state of the cell and transforms it to the next state
            umm{z}=umm{z}-mu.*(network{z}+1)+(mu.*sig)+mu.*input{z};
            u{z}=umm{z}+(beta_e.*input{z});
            cond1=find(network{z}<=0);
            cond2a=(network{z}>0).*(network{z}<(al+u{z})).*(mem{z,(ctime-2)}<=0);
            cond2=find((network{z}>0).*(network{z}<(al+u{z})).*(mem{z,(ctime-2)}<=0));
            cond3=find(((network{z}>=(al+u{z}))+(mem{z,(ctime-2)}>0))>0);
            
            network{z}(cond1)=al./(1-network{z}(cond1))+u{z}(cond1);
            network{z}(cond2)=al+u{z}(cond2); %        memhi{z,ctime}=network{z};
            network{z}(cond3)=-1;
            mem{z,ctime}=network{z};  % records network state
            fired{z}=cond2a; % reset fired to zero
            fcounts{z,cevent}=fcounts{z,cevent}+fired{z};
            fcount{z}=fcount{z}+fired{z};
            tfires=sum(sum(cond2a))+tfires;
            for ffx=1:events
                if ctime==event{ffx}(1)+preset
                    empt=1:(scale{z}(1)*scale{z}(2));
                    if z==1
                        for zz=1:lmax
                            %cellset
                            
                        end
                    end
                end
            end
            if version==2
                input{z}=input{z}*.995;
            elseif version==4
                input{z}=input{z}*.98;
            else
                input{z}=input{z}*.99;
            end
            if version==3 && z==1
                if ctime>1000 && ctime< 2000
                    input{z}=input{z}*0+.05; %#ok<SAGROW>
                end
                if ctime>4000 && ctime< 5000
                    input{z}=input{z}*0+.2; %#ok<SAGROW>
                end
                if ctime>7000 && ctime< 7002
                    input{z}=input{z}*0+.03; %#ok<SAGROW>
                end
            end
            while evec(bean,1)<=ctime
                network{1}(:,evec(bean,2))=network{1}(:,evec(bean,2))+.5;
                bean=bean+1;
            end
            if version<=2
                delaystack{z}=delaystack{z}-countdown;
                for pak=1:stackmax
                    if delaystack{z}(pak,4)==0 && delaystack{z}(pak,5)~=0
                        input{delaystack{z}(pak,1)}(delaystack{z}(pak,2), delaystack{z}(pak,3))=input{delaystack{z}(pak,1)}(delaystack{z}(pak,2),delaystack{z}(pak,3))+delaystack{z}(pak,5);
                    end
                end
            end
        end
        
        for z=1:lmax
            for y=1:scale{z}(1)
                for x=1:scale{z}(2)
                    if fired{z}(y,x)>0;
                        if ctime>50
                            eye((ctime-blink):(ctime+blink))=1;
                        end
                        
                        if stdp==1
                            lastfired{z,y,x}(numel(lastfired{z,y,x})+1)=ctime;
                            
                            
                            for n=1:maxoutputs
                                
                                if connex{z,y,x,n}(1) ~= 0
                                    k=connex{z,y,x,n}(1);
                                    j=connex{z,y,x,n}(2);
                                    i=connex{z,y,x,n}(3);
                                    for nn=1:numel(lastfired{z,y,x})
                                        for oo=1:numel(lastfired{k,j,i})
                                            if abs(lastfired{z,y,x}(nn)-lastfired{k,j,i}(oo))<200
                                                if stdp==1 && connex{z,y,x,n}(5)>0 && lastfired{z,y,x}(nn)<lastfired{k,j,i}(oo) && connex{z,y,x,n}(5)>0
                                                    oconnex=connex{z,y,x,n}(5);
                                                    if k==3
                                                        connex{z,y,x,n}(5)=connex{z,y,x,n}(5)*(1+sscale*exp((lastfired{z,y,x}(nn)-lastfired{k,j,i}(oo))/tau));%innacurate should be changed to 1 to 1 or all to all
                                                        (1+sscale*exp((lastfired{k,j,i}(oo)-lastfired{z,y,x}(nn))/tau));
                                                    end
                                                    if oconnex~=connex{z,y,x,n}(5)
                                                        %connex{z,y,x,n}(5)
                                                        
                                                    end
                                                end
                                                if stdp==1 && connex{z,y,x,n}(5)>0 && lastfired{z,y,x}(nn)>lastfired{k,j,i}(oo) && connex{z,y,x,n}(5)>0
                                                    if k==3
                                                        prev=connex{z,y,x,n}(5);
                                                        connex{z,y,x,n}(5)=connex{z,y,x,n}(5)*(1+wscale*exp((lastfired{k,j,i}(oo)-lastfired{z,y,x}(nn))/tau));
                                                        (1+wscale*exp((lastfired{k,j,i}(oo)-lastfired{z,y,x}(nn))/tau));
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            
                            nn=1;
                            for spike_mem=1:numel(lastfired{z,y,x})
                                if ctime-lastfired{z,y,x}(spike_mem)<200
                                    tlastfired(nn)=lastfired{z,y,x}(spike_mem);
                                    nn=nn+1;
                                end
                            end
                            
                            lastfired{z,y,x}=tlastfired;
                            if n==maxoutputs && j~=0
                                faulty='Dear God, my brain is on fire';
                            end
                        end
                        
                        for n=1:maxoutputs
                            if connex{z,y,x,n}(1)>0
                                k=connex{z,y,x,n}(1);
                                j=connex{z,y,x,n}(2);
                                i=connex{z,y,x,n}(3);
                                s=connex{z,y,x,n}(4);
                                p=connex{z,y,x,n}(5);
                                if n==maxoutputs && p~=0
                                    faulty='Dear God, my brain is on fire';
                                end
                                %%%%%%%%%%%%%%%%%%
                                % STDP
                                %if stdp==1 && connex{z,y,x,n}(5)>0
                                
                                %connex{z,y,x,n}(5)=connex{z,y,x,n}(5)+connex{z,y,x,n}(5)*wscale*exp(lastfired{k}(j,i)-ctime);
                                %end
                                done=0;
                                if version<=2
                                    parker=1;
                                    while done==0
                                        if delaystack{z}(parker,4)<=0
                                            delaystack{z}(parker,:)=connex{z,y,x,n};
                                            done=1;
                                        end
                                        parker=parker+1;
                                    end
                                else
                                    input{k}(j,i)=input{k}(j,i)+p;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %5. graphing
    
    for ctime=1:ttime
        timp1(ctime)=sum(sum(mem{1,ctime}))./numel((mem{1,ctime})); % %#ok<*MSNU>
        timp2(ctime)=sum(sum(mem{2,ctime}))./numel((mem{2,ctime}));
        timp3(ctime)=sum(sum(mem{3,ctime}))/numel((mem{3,ctime}));
        timp4(ctime)=sum(sum(mem{4,ctime}))/numel((mem{4,ctime}));
        %tmpc31(ctime)=mem{3,ctime}(1,1);
        
        
    end
    %     connex{1,1,1,8}
    %     connex{2,1,1,1}
    %     connex{3,1,1,1}
    %     connex{4,1,1,1}
    
    if graph==1
        close all
        chart=figure;
        hold on
        plot(timp1,'r')
        plot(timp2,'b')
        plot(timp3,'g')
        plot(timp4,'y')
        %plot(tmpc31,'m')
        if rem(aeon,1)==0
            saveas(chart,num2str(aeon))
        end
        
    end
    
    responce=ones(3,400)*-.94;
    responcesca=-199:200;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %5. make moviem(plays back network
    if moviem==1 && aeon==aeons
        figure
        for ctime=1:ttime
            if eye(ctime)==1
                if (ctime/1)==int64(ctime/1)
                    for l=1:lmax
                        subplot(2,2,l)
                        arf=image((mem{l,ctime}+1)*200);title(num2str(ctime))
                    end
                    pause(0.02)
                end
            end
        end
    end
    %toc
    
    
    for tunin=1:3
        for flax=1:events
            blitzen=1;
            if flax==events
                capp=ttime;
            else
                capp=event{flax+1}(1);
            end
            for fnix=event{flax}(1):capp
                if version==4
                    trump(blitzen)=mem{4,fnix}((tunin*4-3),1);
                else
                    trump(blitzen)=mem{2,fnix}((tunin*4-3),1);
                end
                blitzen=blitzen+1;
            end
            responce(tunin,(event{flax}(2)+200))=max(trump);
            
        end
        
    end
    figure
    hold on
    plot(responcesca,responce(1,:))
    plot(responcesca,responce(2,:),'r')
    plot(responcesca,responce(3,:),'g')
    rasteraudi
    if stdp==0
        outfoc=4;
        if version==2
            outfoc=2;
        end
        heatchart
        %hgsave(num2str(tilt))
    end
    if version==2
        tiltcalc
    end
    
end

%hotzplotz










