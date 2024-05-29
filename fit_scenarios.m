clear all
close all

ss = @(x) x/x(5);


exportin = [1.246, 1.112, 1.013,1.057, 1.0,1.0];
Eflat = 0.7326;

areas = [12.56	11.22	10.82	10.33	8.794 9.539	];

u0vals = [1.05, .925, .85, .8, .72, .75];



kvals = [2.11, 1.747, 1.412, 1.56, 1];

Lvals = [ 0.3640	0.4973	0.5446	0.7003	0.7089, 0.7089];


flux_data = [2.8, 1.9, 1.44, 1.039, 1.0, .9601];
flux_data_errors = [.175, .114, .067, .08466, .07237, .08178];


simple_model_flux = @(D, k, L, u0, E) E.*D.*k.*u0./(D+k.*L);

simple_model_flux_s = @(D, k, L, u0, E) ss(E.*D.*k.*u0./(D+k.*L));


loss_u0exp = @(lT) mean((ss(simple_model_flux(exp(lT(1)),exp(lT(2)),Lvals(1:5), u0vals(1:5), exportin(1:5)))-flux_data(1:5)).^2);


LB = [-2,-2];
UB = [2,2];
x0=[0,0];

[x1,fval1,exitflag1,output1] = fminsearchbnd(loss_u0exp,x0,LB,UB);



%u0colors  =  multigradient([100, 143, 255;255,176,0]/256,length=4);


figure;


%plot(1:5, simple_model_flux_s(exp(x1(1)), exp(x1(2)), Lvals, u0vals, exportin ),'DisplayName','full model, u0 decreasing ');
%hold on;
%plot(1:5, simple_model_flux_s(exp(x2(1)), exp(x2(2)), Lvals, 1, exportin ),'DisplayName','full model, u0 const');
%plot(1:5, simple_model_flux_s(exp(x3(1)), exp(x3(2)), Lvals, u0flipped1, exportin ),'DisplayName','full model, u0 increas');
%plot(1:5, simple_model_flux_s(exp(x4(1)), exp(x4(2)), Lvals, u0flipped2, exportin ),'DisplayName','u0 increasing more');



 load('sweep_geometries15-Feb-2024-22-18-45.mat')

 nondim_scale = 10;

k =.7324*nondim_scale; % from best fit file

Dfast=1*nondim_scale; 
Dslow =.0498*nondim_scale; % from best fit file

decay=0;
box_width=15;
box_height=5;
hmax=.1;


fluxes_out_full = zeros(6,1);

parfor i =  1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    kval = k;
    kflat = k;
    x_nuc = squeeze(xvals_chrom(i,:,:));
    

  e_i = exportin(i);
    
    flux_i_full = solve_pde_flux_out(box_width,box_height,pillar_rad,pillar_dist,...
        pillar_height,...
        x_nuc,u0,kval,kflat, decay,Dfast,Dslow, cthick,cthick_nuc,e_i, Eflat,hmax,0);  

    fluxes_out_full(i) = flux_i_full;

close all;
end 

colors = [100, 143, 255; ...
    120,94 ,240;...
    220 ,38 ,127;...
    254 ,97 ,0;...
    255 ,176 ,0]/256;

flux_simp = simple_model_flux(exp(x1(1)), exp(x1(2)), Lvals, u0vals, exportin );

p_indep = simple_model_flux(exp(x1(1)), exp(x1(2)), Lvals, u0vals(5), exportin )/flux_simp(5);
L_indep = simple_model_flux(exp(x1(1)), exp(x1(2)), Lvals(5), u0vals, exportin )/flux_simp(5);
E_indep = simple_model_flux(exp(x1(1)), exp(x1(2)), Lvals, u0vals, exportin(5) )/flux_simp(5);
PL_indep = simple_model_flux(exp(x1(1)), exp(x1(2)), Lvals(5), u0vals(5), exportin )/flux_simp(5);

%hold on;
hold on;

plot(fliplr(p_indep),'diamond-','MarkerSize',10,'LineWidth',1.5,'Color',colors(2,:),'displayname','p independent');
plot(fliplr(L_indep),'^-','MarkerSize',10,'LineWidth',1.5,'Color',colors(3,:),'displayname','L independent');
plot(fliplr(E_indep),'<-','MarkerSize',10,'LineWidth',1.5,'Color',colors(4,:),'displayname','E independent');
plot(fliplr(PL_indep),'v-','MarkerSize',10,'LineWidth',1.5,'Color',colors(5,:),'displayname','p+L independent');
plot(fliplr(ss(fluxes_out_full)'),'o--','MarkerSize',10,'LineWidth',1.5,'Color',colors(1,:),'displayname','full model')
plot(fliplr(simple_model_flux_s(exp(x1(1)), exp(x1(2)), Lvals, u0vals, exportin )),'square-','MarkerSize',10,'LineWidth',1.5,'Color',colors(1,:),'displayname','approx model')



legend


 errorbar(1:6,[fliplr(flux_data)], [fliplr(flux_data_errors)],'o','MarkerSize', 8,'MarkerEdgeColor',[25,25,25]/256,'DisplayName','RPL13 data','LineWidth',1.5,...
     'Color','black')
 box off

 xticks([1,2,3,4,5,6,7])
 labels = fliplr({"R150","R200","R300","R300","R350","R400"});
 xticklabels(labels)
 ylabel('normalized flux out')
 xlim([0.5 6.5])
 legend box off;
 set(gca,'FontSize',14)
set(gca,'LineWidth',1.25)


%figure;

%scatter(1:6,flipud(fluxes_out_full),50,fluxes_out_full,'filled')
%colorbar
%colormap('cmasher')