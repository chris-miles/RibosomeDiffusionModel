clear all
close all

 load('sweep_geometries15-Feb-2024-22-18-45.mat')

fluxes_out_full = zeros(6,1);
fluxes_simple= zeros(6,1);


flux_data = [ 2.8, 1.9, 1.44, 1.039, 1.0, NaN];
flux_data_errors = [ .175, .114, .067, .08466, .07237,NaN];

exportin = [ 1.246, 1.112, 1.013,1.057, 1.0, 1.0];
Eflat = 0.7326;


u0vals = [ 1.05, .925, .85, .8, .72, .75];



nondim_scale = 1;

k= 0.7324*nondim_scale;
kflat = k;


Dfast=1*nondim_scale;
Dslow =.0498*nondim_scale;

decay=0;
box_width=15;
box_height=5;
hmax=.25;

Lvals = [0.3640	0.4973	0.5446	0.7003	0.7089 0.7089];

simple_model_flux = @(D, k, L, u0) D*k*u0/(D+k*L);

for i =  1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    kval = k;
    x_nuc = squeeze(xvals_chrom(i,:,:));
    e_i = exportin(i);
    

    flux_i_full = solve_pde_flux_out(box_width,box_height,pillar_rad,pillar_dist,...
        pillar_height,...
        x_nuc,u0,kval,kflat, decay,Dfast,Dslow, cthick,cthick_nuc,e_i, Eflat,hmax,1);  

 
    fluxes_out_full(i) = flux_i_full;

   
end


fluxes_out_full = fluxes_out_full/fluxes_out_full(end-1);




 figure;
 plot(fluxes_out_full,'o-','LineWidth',2,'Color',[100, 143, 255]/256,'DisplayName','full model',...
     'MarkerSize',8)
 hold on;



 errorbar(1:6,flux_data(1:6), flux_data_errors(1:6),'o','MarkerSize', 8,'MarkerEdgeColor',[25,25,25]/256,'DisplayName','RPL13 data','LineWidth',1.5,...
     'Color','black')

 xticks([1,2,3,4,5,6,7])
 labels = {"R150","R200","R300","R300","R350","R400"};
 xticklabels(labels)
 ylabel('normalized flux out')
 
 legend box off;
 set(gca,'FontSize',14)
set(gca,'LineWidth',1.25)
