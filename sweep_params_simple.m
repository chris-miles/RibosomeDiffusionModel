clear all
close all


fluxes_Dslow_sweep = zeros(6,11);
fluxes_k_sweep = zeros(6,11);

kmeasure = [ 2.11, 1.747, 1.412, 1.56, 1, 1];

Lvals = [0.3640	0.4973	0.5446	0.7003	0.7089 0.7089];

exportin = [1.246, 1.112, 1.013,1.057, 1.0, 1.0];

u0vals = [1.05, .925, .85, .8, .72, .75];

nondim_scale=10;

ksimp=.2691*nondim_scale;
Dslowsimp =.0698*nondim_scale;

simple_model_flux = @(D, k, L, u0,E) D*k*u0/(D+k*L);

%colors  =  multigradient([5,113,176;155,155,155;202,0,32]/256,length=11);
colors  =  multigradient([0,90,181;155,155,155;220,50,32]/256,length=11);


Dslowvals = logspace(-2,2,11)*Dslowsimp;
kvals = logspace(-2,2,10)*ksimp;

for i = 1:6

    u0_i = u0vals(i);
    cthick_i = Lvals(i);
    exp_i = exportin(i);
    
    for j = 1:10
flux_ij_Dslow = simple_model_flux(Dslowvals(j),ksimp, cthick_i, u0_i, exp_i);
flux_ij_k = simple_model_flux(Dslowsimp,kvals(j), cthick_i, u0_i, exp_i);

fluxes_Dslow_sweep(i,j) = flux_ij_Dslow;
fluxes_k_sweep(i,j) = flux_ij_k;
    end 
end 



figure('position',[0,0,900,700]);
subplot(2,2,1);
hold on;
for i = 1:10
plot(fluxes_Dslow_sweep(:,i),'o-','color',colors(i,:),'LineWidth',2)
end 
set(gca,'ColorScale','log')
labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 7])
 xticklabels(labels)
  xlim([0.5 6.5])
 %axis square;
 set(gca,'YScale','log')
 set(gca,'LineWidth',1.5)
 ylabel('unnormalized flux')
 title('D_{slow} varied','FontWeight','normal')
 ylim([1e-3 1e1])


subplot(2,2,3);
hold on;
for i = 1:10
plot(fluxes_Dslow_sweep(:,i)/fluxes_Dslow_sweep(end-1,i),'o-','color',colors(i,:),'LineWidth',2)
end 
set(gca,'ColorScale','log')
labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 7])
 xticklabels(labels)
  xlim([0.5 6.5])
 %axis square;
 set(gca,'LineWidth',1.5)
 ylabel('flux normalized by R350')
 title('D_{slow} varied','FontWeight','normal')
 ylim([0.9 3])


subplot(2,2,2);
hold on;
for i = 1:10
plot(fluxes_k_sweep(:,i),'o-','color',colors(i,:),'LineWidth',2)
end 
set(gca,'ColorScale','log')
labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 7])
 xticklabels(labels)
  xlim([0.5 6.5])
 %axis square;
  set(gca,'YScale','log')
 set(gca,'LineWidth',1.5)
 ylabel('unnormalized flux')
 title('k varied','FontWeight','normal')
 ylim([1e-3 10^1])

subplot(2,2,4);
hold on;
for i = 1:10
plot(fluxes_k_sweep(:,i)/fluxes_k_sweep(end-1,i),'o-','color',colors(i,:),'LineWidth',2)
end 
set(gca,'ColorScale','log')
labels = {"R150","R200","R300","R300","R350","R400"};
 xticks([1 2 3 4 5 6 7])
 xticklabels(labels)
  xlim([0.5 6.5])
 %axis square;
 set(gca,'LineWidth',1.5)
 ylabel('flux normalized by R350')
 title('k varied','FontWeight','normal')
 ylim([0.9 3])

 