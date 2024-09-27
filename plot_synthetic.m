clear all,
close all,
format long,
dt=2.5/8191;
offsets=[130 330 520 670 770 870 1020 1100 1230 1330 1450 1600 1750 1850 2000 2100 2250 2350 2450];
disp('swich=1: Plot of P-SV transmission coeffecient with incident angle');
disp('swich=2: Plot of P-SV transmission coeffecient with incident angle');
disp('swich=3:  Plot incident angle versus depth for each angle');
disp('swich=4: save time and amplitudes of each trace in an individual txt file');
disp('swich=5: Calculate synthetic time of P and S waves');
disp('swich=6: Plot synthetic seismograms');
swich = input('swich=');
 vp1=5000; ro1=3; vs1=vp1/sqrt(3); % this is properties of layer1
    vp2=3265.2; ro2=2.335; vs2=vp2/sqrt(3); % this is properties of layer2
    i=0:0.01:pi/2;
    p=sin(i)./vp1;
    q=2*(ro2*(vs2^2)-ro1*(vs1^2));
    X=ro2-q*(p.^2);
    Y=ro1+q*(p.^2);
    Z=ro2-ro1-q*(p.^2);
    P1=((1-(vp1.^2)*(p.^2)).^0.5);
    P2=((1-(vs1.^2)*(p.^2)).^0.5);
    P3=((1-(vp2.^2)*(p.^2)).^0.5);
    P4=((1-(vs2.^2)*(p.^2)).^0.5);
    D=((q^2).*(p.^2).*P1.*P2.*P3.*P4)+(ro1*ro2.*(vs1*vp2.*P1.*P4+vp1*vs2.*P2.*P3))+(vp1*vs1.*P3.*P4.*(Y.^2))+(vp2*vs2.*P1.*P2.*(X.^2))+(vp1*vp2*vs1*vs2*(p.^2).*(Z.^2));% P-P reflection
    Dinv=(1./D);
%==============================================================================================
%=================  Plot of P-SV reflection1 coeffecient with incident angle ================= 
if (swich==1)
    R33 = Dinv.*(q^2*p.^2.*P1.*P2.*P3.*P4 + ro1*ro2*(vs1*vp2*P1.*P4-vp1*vs2*P2.*P3)-vp1*vs1*P3.*P4.*Y.^2 + vp2*vs2*P1.*P2.*X.^2-vp1*vp2*vs1*vs2*p.^2.*Z.^2);
    R31=2*vp1.*p.*P1.*Dinv.*((q.*P3.*P4.*Y)+(vp2*vs2.*X.*Z)); % PSV reflection
i=rad2deg(i);
figure,
plot(i,R31,'r-','linewidth',2);
xlabel('Incident Angle [degree]','fontsize',35);ylabel('Reflection Coeff','fontsize',35);
set(gca,'fontsize',30,'Linewidth',3,'Ytick',0:0.25:3);box on;title('Reflection P-SV','fontsize',40);
end
%==============================================================================================
%==============  Plot of P-SV transmission coeffecient with incident angle ====================
if (swich==2)
T31=2*vp1*ro1*p.*P1.*Dinv.*(q.*P2.*P3-vp2*vs1*Z); % Transmission  coeff  P-SV 
i=rad2deg(i);
figure,
plot(i,T31,'r-','linewidth',2);
xlabel('Incident Angle [degree]','fontsize',35);ylabel('Transmission Coeff','fontsize',35);
set(gca,'fontsize',30,'Linewidth',3,'Ytick',0:0.25:3);box on;title('Transmission P-SV','fontsize',40);
end
%==============================================================================================
%==============================  Plot incident angle versus depth  ============================
if (swich==3)
depths=(1145:20:1145+20*69)'
offsets=offset:100:offset;
for i=1:1%length(offsets)
figure, 
hold on,
plot(depths,rad2deg(atan(offsets(i)./depths)),'ko');title(['Offsets=',num2str(offsets(i)),' m'],'fontsize',40);
plot([min(depths) max(depths)], [10 10],'r-','linewidth',3);
plot([min(depths) max(depths)], [25 25],'r-','linewidth',3);
ylabel('Incident Angle [degree]','fontsize',35);xlabel('Depth [m]','fontsize',35);
set(gca,'fontsize',30,'Linewidth',3);box on;
end
end
nt=8192;      %%%% sampling rate %%%%%%
%===============================================================================================
%============= save time and amplitudes of each trace in an individual txt file ================
if (swich==4)
for k=1:1%length(offsets)
amplx=load(['Synthetic_atten\offset',num2str(offsets(k)),'\deplacx']);
amply=load(['Synthetic_atten\offset',num2str(offsets(k)),'\deplacy']);
amplz=load(['Synthetic_atten\offset',num2str(offsets(k)),'\deplacz']);
% mkdir(['Synthetic_atten\offset',num2str(offsets(k)),'\Horizx_source']);
% mkdir(['Synthetic_atten\offset',num2str(offsets(k)),'\Horizy_source']);
% mkdir(['Synthetic_atten\offset',num2str(offsets(k)),'\vertical_source']);
   for i=nt:nt:length(amplx)
      vx=amplx(i-nt+1:i,:);vy=amply(i-nt+1:i,:);vz=amplz(i-nt+1:i,:);
      dlmwrite(['Synthetic_atten\offset',num2str(offsets(k)),'\Horizx_source\tracex',num2str(i/nt),'.txt'],vx,'delimiter',' ');
      dlmwrite(['Synthetic_atten\offset',num2str(offsets(k)),'\Horizy_source\tracey',num2str(i/nt),'.txt'],vy,'delimiter',' ');
      dlmwrite(['Synthetic_atten\offset',num2str(offsets(k)),'\vertical_source\tracez',num2str(i/nt),'.txt'],vz,'delimiter',' ');
   end
end
end
%===============================================================================================
%============================ Calculate synthetic time of P and S waves ========================
if (swich==5)
 input=load(['Synthetic_atten\offset',num2str(offsets(1)),'\skb_mod.txt']);
 depth=[0;input(:,1)*1000];vp=input(:,2)*1000;vs=input(:,3)*1000;
% depth(2)=0;
 for  i=1:length(offsets)
 fid=fopen(['Synthetic_atten\synt_time_offset',num2str(offsets(i))],'w');
 for k=5:length(depth)-2   % P waves and PS converted at sea floor
     theta=atan(offsets(i)/depth(k));
     tp=sum(diff(depth(1:4))./(cos(theta)*vp(1:3)))+sum(diff(depth(4:k))./(cos(theta)*vp(4:k-1)));
     tps1=sum(diff(depth(1:2))./(cos(theta)*vp(1)))+sum(diff(depth(2:4))./(cos(theta)*vs(2:3)))+sum(diff(depth(4:k))./(cos(theta)*vs(4:k-1)));
     tps2=0;
      if (k>52)  % PS Coverted at Um_Nahr layer
      tps0=sum(diff(depth(1:4))./(cos(theta)*vp(1:3)))+sum(diff(depth(4:52))./(cos(theta)*vp(4:52-1)));
      tps2=tps0+sum(diff(depth(52:k))./(cos(theta)*vs(52:k-1))); 
      end
     fprintf(fid,'\t%f\t%f\t%f\t%f\t%f \n',depth(k),tp,tps1,tps2,rad2deg(theta));
 end
%  for k=5:length(depth)-2
%  theta=atan(offset/depth(k));
%  t_upp_PP=sum(diff(depth(1:4))./(cos(theta)*vp(1:3)));t_upp_PS=diff(depth(1:2))/(cos(theta)*vs(1))+sum(diff(depth(2:4))./(cos(theta)*vs(2:3)));t_upp_SS=sum(diff(depth(1:4))./(cos(theta)*vs(1:3)));
%  tppp=diff(depth(4:k))./vp(4:k-1);tpp=t_upp_PP+sum(tppp)/cos(theta);
%  tpps=diff(depth(4:k))./vs(4:k-1);tps=t_upp_PS+sum(tpps)/cos(theta);
%  tsss=diff(depth(4:k))./vs(4:k-1);tss=t_upp_SS+sum(tsss)/cos(theta);
%  fprintf(fid,'\t%f\t%f\t%f\t%f \n',depth(k),tpp+t0,tps+t0,tss+t0);
%  end
fclose(fid);
 end
end
%===============================================================================================
%================================= Plot synthetic seismograms  ==================================
if (swich==6)
 clear i
 offsets=[130 330 520 670 770 870 1020 1100 1230 1330 1450 1600 1750 1850 2000 2100 2250 2350 2450];
depth=1145:20:2365;
n=62;nt=8192;t0=0.08;index=21;
trac1=load(['Synthetic2w_atten\offset',num2str(offsets(1)),'\Vertical_source\tracez1.txt']); % Number of traces
t=trac1(:,1);tmin=0; tmax=1.5; % minimum and maximum time
for i=1:1%length(offsets)
syn_time=load(['Synthetic2w_atten\synt_time_offset',num2str(offsets(i))]);
seismz=load(['estimate_Q\refl\seismo',num2str(offsets(i)),'_Hz.txt']);
%x=depth(index):20:depth(end);
x=syn_time(index:end,5);y=t-0.08;
figure,
hold on,
wiggle(y',x,seismz(index:end,:)')
hold on,
plot(syn_time(index:end,5),syn_time(index:end,2),'r-','Linewidth',1);
plot(syn_time(index:end,5),syn_time(index:end,3),'g-','Linewidth',1);
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',round(min(x)):round((max(x)-min(x))/7):round(max(x)),...
    'XTicklabel',round(min(x)):round((max(x)-min(x))/7):round(max(x)),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);
box on;axis([min(x)-0.25 max(x)+0.25 0 2]);
xlabel('Depth [m]','fontsize',11);ylabel('Time [sec]','fontsize',11);
    title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   
end
end
%============================================================================================================================
%================================= Plot synthetic seismograms and save them in a pppt file ==================================
if (swich==7)
 clear i
 offsets=[130 330 520 670 770 870 1020 1100 1230 1330 1450 1600 1750 1850 2000 2100 2250 2350 2450];
% Add image in a vector (non-raster) format
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
exportToPPTX('new','Dimensions',[12 6], ...
    'Title','Example Presentation', ...
    'Author','MatLab', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');

depth=1145:20:2365;
n=62;nt=8192;index=1;
trac1=load(['Synthetic_atten\offset',num2str(offsets(1)),'\Vertical_source\tracez1.txt']); % Number of traces
t=trac1(:,1);tmin=0; tmax=1.5;t0=0.08; % minimum and maximum time
x=depth(index):20:depth(49);y=t-t0;
trac1=load(['Synthetic_atten\offset',num2str(offsets(1)),'\Vertical_source\tracez1.txt']); % Number of traces
t=trac1(:,1);tmin=0; tmax=1.5;depthmin=1145+(index-1)*20;depthmax=2365; % minimum and maximum time
for i=1:length(offsets)
syn_time=load(['Synthetic_atten\synt_time_offset',num2str(offsets(i))]);
% seismx=load(['estimate_Q\refl\seismo',num2str(offsets(i)),'_Hx.txt']);
% seismy=load(['estimate_Q\refl\seismo',num2str(offsets(i)),'_Hy.txt']);
% seismz=load(['estimate_Q\refl\seismo',num2str(offsets(i)),'_Hz.txt']);
seismxd=load(['estimate_Q\No_refl\seismo',num2str(offsets(i)),'_Hx.txt']);
seismyd=load(['estimate_Q\No_refl\seismo',num2str(offsets(i)),'_Hy.txt']);
seismzd=load(['estimate_Q\No_refl\seismo',num2str(offsets(i)),'_Hz.txt']);

% figure,
% hold on,
% wiggle(y',x,seismx(index:end,:)')
% hold on,
% plot(syn_time(index:end,1),syn_time(index:end,2),'r-');
% plot(syn_time(index:end,1),syn_time(index:end,3),'g-');
% plot(syn_time(53:end,1),syn_time(53:end,4),'b-');
% set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',depth(index):200:depth(end),...
%     'XTicklabel',depth(index):200:depth(end),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);box on;
% axis([depth(index) depth(end)+20 0 2]);
% xlabel('Depth [m]','fontsize',11);ylabel('Time [sec]','fontsize',11);title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   
% 
%     saveas(gcf,'vectorFile','png');
% slideNum = exportToPPTX('addslide');
% exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
% exportToPPTX('addtext',sprintf('Downgoin+upgoing: Hx component'),'Position',[0.54 2.42 2.29 0.7]);
% close(gcf); 
% 
% delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismxd(index:49,:)')
hold on,
plot(syn_time(index:49,1),syn_time(index:49,2),'r-');
plot(syn_time(index:49,1),syn_time(index:49,3),'g-');
% plot(syn_time(53:end,1),syn_time(53:end,4),'b-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',depth(index):200:depth(49),...
    'XTicklabel',depth(index):200:depth(49),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);box on;
axis([depth(index) depth(49)+20 0 2]);
xlabel('Depth [m]','fontsize',11);ylabel('Time [sec]','fontsize',11);title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoin waves: Hx component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 

delete('vectorFile.png');

% figure,
% hold on,
% wiggle(y',x,seismy(index:end,:)')
% hold on,
% plot(syn_time(index:end,1),syn_time(index:end,2),'r-');
% plot(syn_time(index:end,1),syn_time(index:end,3),'g-');
% plot(syn_time(53:end,1),syn_time(53:end,4),'b-');
% set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',depth(index):200:depth(end),...
%     'XTicklabel',depth(index):200:depth(end),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);box on;
% axis([depth(index) depth(end)+20 0 2]);
% xlabel('Depth [m]','fontsize',11);ylabel('Time [sec]','fontsize',11);title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   
% 
%     saveas(gcf,'vectorFile','png');
% slideNum = exportToPPTX('addslide');
% exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
% exportToPPTX('addtext',sprintf('Downgoing+upgoing: Hy component'),'Position',[0.54 2.42 2.29 0.7]);
% close(gcf); 
% 
% delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismyd(index:49,:)')
hold on,
plot(syn_time(index:49,1),syn_time(index:49,2),'r-');
plot(syn_time(index:49,1),syn_time(index:49,3),'g-');
plot(syn_time(53:49,1),syn_time(53:49,4),'b-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',depth(index):200:depth(49),...
    'XTicklabel',depth(index):200:depth(49),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);box on;
axis([depth(index) depth(49)+20 0 2]);
xlabel('Depth [m]','fontsize',11);ylabel('Time [sec]','fontsize',11);title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoin waves: Hy component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 

delete('vectorFile.png');

% figure,
% hold on,
% wiggle(y',x,seismz(index:end,:)')
% hold on,
% plot(syn_time(index:end,1),syn_time(index:end,2),'r-');
% plot(syn_time(index:end,1),syn_time(index:end,3),'g-');
% plot(syn_time(53:end,1),syn_time(53:end,4),'b-');
% set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',depth(index):200:depth(end),...
%     'XTicklabel',depth(index):200:depth(end),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);box on;
% axis([depth(index) depth(end)+20 0 2]);
% xlabel('Depth [m]','fontsize',11);ylabel('Time [sec]','fontsize',11);title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   
% 
%     saveas(gcf,'vectorFile','png');
% slideNum = exportToPPTX('addslide');
% exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
% exportToPPTX('addtext',sprintf('Downgoing+upgoing: Hz component'),'Position',[0.54 2.42 2.29 0.7]);
% close(gcf); 
% 
% delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismzd(index:49,:)')
hold on,
plot(syn_time(index:49,1),syn_time(index:49,2),'r-');
plot(syn_time(index:49,1),syn_time(index:49,3),'g-');
% plot(syn_time(53:end,1),syn_time(53:end,4),'b-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',depth(index):200:depth(49),...
    'XTicklabel',depth(index):200:depth(49),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);box on;
axis([depth(index) depth(49)+20 0 2]);
xlabel('Depth [m]','fontsize',11);ylabel('Time [sec]','fontsize',11);title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoin waves: Hz component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 
delete('vectorFile.png');
end

% Filename automatically checked for proper extension 
newFile = exportToPPTX('saveandclose','synthetic_seismo_conv_NU_Dow'); 
end
%===============================================================================================
%================================= Plot synthetic seismograms versus incident angle and save them in a ppt file ==================================
if (swich==8)
 clear i
 offsets=[130 330 520 670 770 870 1020 1100 1230 1330 1450 1600 1750 1850 2000 2100 2250 2350 2450];
% Add image in a vector (non-raster) format
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
exportToPPTX('new','Dimensions',[12 6], ...
    'Title','Example Presentation', ...
    'Author','MatLab', ...
    'Subject','Automatically generated PPTX file', ...
    'Comments','This file has been automatically generated by exportToPPTX');

n=62;nt=8192;index=21;
trac1=load(['Synthetic2w_atten\offset',num2str(offsets(1)),'\Vertical_source\tracez1.txt']); % Number of traces
t=trac1(:,1);y=t-0.08;tmin=0; tmax=1.5;t0=0.08; % minimum and maximum time
%x=depth(index):20:depth(end);y=t-t0;
trac1=load(['Synthetic2w_atten\offset',num2str(offsets(1)),'\Vertical_source\tracez1.txt']); % Number of traces
t=trac1(:,1);tmin=0; tmax=1.5;depthmin=1145+(index-1)*20;depthmax=2365; % minimum and maximum time
for i=1:length(offsets)
syn_time=load(['Synthetic2w_atten\synt_time_offset',num2str(offsets(i))]);
x=syn_time(index:end,5);
seismx=load(['estimate_Q\refl\seismo',num2str(offsets(i)),'_Hx.txt']);
seismy=load(['estimate_Q\refl\seismo',num2str(offsets(i)),'_Hy.txt']);
seismz=load(['estimate_Q\refl\seismo',num2str(offsets(i)),'_Hz.txt']);
seismxd=load(['estimate_Q\No_refl\seismo',num2str(offsets(i)),'_Hx.txt']);
seismyd=load(['estimate_Q\No_refl\seismo',num2str(offsets(i)),'_Hy.txt']);
seismzd=load(['estimate_Q\No_refl\seismo',num2str(offsets(i)),'_Hz.txt']);

figure,
hold on,
wiggle(y',x,seismx(index:end,:)')
hold on,
plot(syn_time(index:end,5),syn_time(index:end,2),'r-');
plot(syn_time(index:end,5),syn_time(index:end,3),'g-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),...
    'XTicklabel',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);
box on;axis([min(x)-abs(x(1)-x(2)) max(x)+abs(x(1)-x(2)) 0 2]);
xlabel('Incid-angl [^o]','fontsize',11);ylabel('Time [sec]','fontsize',11);
    title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoin+upgoing: Hx component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 

delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismxd(index:end,:)')
hold on,
plot(syn_time(index:end,5),syn_time(index:end,2),'r-','Linewidth',1);
plot(syn_time(index:end,5),syn_time(index:end,3),'g-','Linewidth',1);
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),...
    'XTicklabel',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);
box on;axis([min(x)-abs(x(1)-x(2)) max(x)+abs(x(1)-x(2)) 0 2]);
xlabel('Incid-angl [^o]','fontsize',11);ylabel('Time [sec]','fontsize',11);
    title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);    

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoin waves: Hx component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 

delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismy(index:end,:)')
hold on,
plot(syn_time(index:end,5),syn_time(index:end,2),'r-');
plot(syn_time(index:end,5),syn_time(index:end,3),'g-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),...
    'XTicklabel',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);
box on;axis([min(x)-abs(x(1)-x(2)) max(x)+abs(x(1)-x(2)) 0 2]);
xlabel('Incid-angl [^o]','fontsize',11);ylabel('Time [sec]','fontsize',11);
    title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);    

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoing+upgoing: Hy component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 

delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismyd(index:end,:)')
hold on,
plot(syn_time(index:end,5),syn_time(index:end,2),'r-');
plot(syn_time(index:end,5),syn_time(index:end,3),'g-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),...
    'XTicklabel',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);
box on;axis([min(x)-abs(x(1)-x(2)) max(x)+abs(x(1)-x(2)) 0 2]);
xlabel('Incid-angl [^o]','fontsize',11);ylabel('Time [sec]','fontsize',11);
    title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoin waves: Hy component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 

delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismz(index:end,:)')
hold on,
plot(syn_time(index:end,5),syn_time(index:end,2),'r-');
plot(syn_time(index:end,5),syn_time(index:end,3),'g-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),...
    'XTicklabel',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);
box on;axis([min(x)-abs(x(1)-x(2)) max(x)+abs(x(1)-x(2)) 0 2]);
xlabel('Incid-angl [^o]','fontsize',11);ylabel('Time [sec]','fontsize',11);
    title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);     

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoing+upgoing: Hz component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 

delete('vectorFile.png');

figure,
hold on,
wiggle(y',x,seismzd(index:end,:)')
hold on,
plot(syn_time(index:end,5),syn_time(index:end,2),'r-');
plot(syn_time(index:end,5),syn_time(index:end,3),'g-');
set(gca,'fontsize',10,'Linewidth',2,'XMinorTick','on','YMinorTick','on','XTick',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),...
    'XTicklabel',round(min(x),2):round((max(x)-min(x))/7,2):round(max(x),2),'YTick',tmin:0.3:tmax,'YTicklabel',tmin:0.3:tmax);
box on;axis([min(x)-abs(x(1)-x(2)) max(x)+abs(x(1)-x(2)) 0 2]);
xlabel('Incid-angl [^o]','fontsize',11);ylabel('Time [sec]','fontsize',11);
    title(['Offset=',num2str(offsets(i)),'m'],'fontsize',16);   

    saveas(gcf,'vectorFile','png');
slideNum = exportToPPTX('addslide');
exportToPPTX('addpicture','vectorFile.png','Position',[2.83 0 8 6])
exportToPPTX('addtext',sprintf('Downgoin waves: Hz component'),'Position',[0.54 2.42 2.29 0.7]);
close(gcf); 
delete('vectorFile.png');
end

% Filename automatically checked for proper extension 
newFile = exportToPPTX('saveandclose','synthetic_seismo_inc_angles'); 
end