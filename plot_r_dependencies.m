clear all
close all


exportin = [1.246, 1.112, 1.013,1.057, 1.0, 1.0];

areas = [12.56	11.22	10.82	10.33	8.794	9.539	];

u0vals = [1.05, .925, .85, .8, .72, .75];

u0flipped1 = linspace(.72,1.05,5);
u0flipped2 = linspace(.72/2,1.05,5);


kflat = 1.138;
Eflat = .7326;

kvals = [ 2.11, 1.747, 1.412, 1.56, 1, 1];

Lvals = [ 0.3640	0.4973	0.5446	0.7003	0.7089];
Lflat = .73;

figure('position',[500,500,1000,650]);


u0colors  =  multigradient([100, 143, 255;255,176,0]/256,length=4);


subplot(2,3,1);
plot(u0vals,'o-','MarkerSize',10,'LineWidth',2,'Color',u0colors(1,:),'displayname','p decreasing (measured)')
hold on;

plot(mean(u0vals)*ones(1,6),'diamond-','MarkerSize',10,'LineWidth',2,'Color',u0colors(2,:),'displayname','p independent');
%plot(u0flipped1,'^-','MarkerSize',10,'LineWidth',2,'Color',u0colors(3,:),'displayname','p increasing');
%plot(u0flipped2,'square-','MarkerSize',10,'LineWidth',2,'Color',u0colors(4,:),'displayname','p increasing+');

 labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 ])
 xticklabels(labels)
 title({'nuclear pre-ribosome', 'production, p'},'FontWeight','normal')
 set(gca,'LineWidth',1.5)
 box off;
 axis square;
xlim([0 7])
legend('location',"best")
legend box off;

subplot(2,3,2);
plot( Lvals,'o-','MarkerSize',10,'LineWidth',2,'Color',[100, 143, 255]/256,'DisplayName','measured')
hold on;
yline(Lflat,'DisplayName','L_{flat} (measured)','LineWidth',2,'LineStyle','--','Color',[0.5,0.5,0.5])

 labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 ])
 xticklabels(labels)
 title({'heterochromatin thickness','on pillars, L [um]'},'FontWeight','normal')
 set(gca,'LineWidth',1.5)
 box off;
 axis square;
 legend('location',"best")
legend box off;

xlim([0 7])
ylim([0.3 0.8])

subplot(2,3,3);
plot( kvals,'x-','MarkerSize',10,'LineWidth',2,'Color',[254, 97, 0]/256, 'displayname','k increasing (measured)')
hold on;
plot(mean(kvals)*ones(1,6),'o-','MarkerSize',10,'LineWidth',2,'Color',u0colors(1,:), 'displayname','NPC saturated');
yline(kflat,'DisplayName','k_{flat} (measured)','LineWidth',2,'LineStyle','--','Color',[0.5,0.5,0.5])

 labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 ])
 xticklabels(labels)
 title({'surface reactivity','~ NPC concentration, k'},'FontWeight','normal')
 set(gca,'LineWidth',1.5)
 box off;
 axis square;

xlim([0 7])
ylim([0.9 2.3])
legend('location',"best")
legend box off;

 subplot(2,3,4);
 plot( exportin,'o-','MarkerSize',10,'LineWidth',2,'Color',u0colors(1,:),'displayname','E increasing (measured)')
 hold on;
 plot(mean(exportin)*ones(1,6),'x-','MarkerSize',10,'LineWidth',2,'Color',[254, 97, 0]/256,'displayname','exportin saturated')
yline(Eflat,'DisplayName','E_{flat} (measured)','LineWidth',2,'LineStyle','--','Color',[0.5,0.5,0.5])


 labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 ])
 xticklabels(labels)
 title('exportin, E','FontWeight','normal')
 set(gca,'LineWidth',1.5)
  box off;
  axis square;
xlim([0 7])
ylim([0.65 1.3])
legend('location',"best")
legend box off;

 subplot(2,3,5);
 plot(areas,'o-','MarkerSize',10,'LineWidth',2,'Color',[100, 143, 255]/256)
 labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 ])
 xticklabels(labels)
 title('nucleolar areas [um^2]','FontWeight','normal')
 set(gca,'LineWidth',1.5)
 box off;
 axis square; 
 xlim([0 7])
 ylim([ 8 13])