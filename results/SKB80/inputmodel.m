%         ========   PLoet seismograms ========
uz=load('deplacz');
nb_samp=1024;depth_rec=2330:8.8:2550;
uz_resh=reshape(uz(:,2),[nb_samp 26]);
figure,
wiggle((1:1024)*(1/1024)',depth_rec',uz_resh)
set(gca,'LineWidth',2,'XTick',2300:10:2600,'XTickLabel',2300:20:2600,'YTick',0:0.2:1,'YTickLabel',{(0:0.2:1)},'TickDir','out','TickLength',[0.01 0.01],'XMinorTick','on','YMinorTick','on','Fontsize',25);box on;
 xlabel('depth [m]','Fontsize',35);ylabel('Time [s]','Fontsize',35);axis([2300 2580 0 1]);
