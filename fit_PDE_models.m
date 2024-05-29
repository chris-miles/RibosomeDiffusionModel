clear all;
close all;


load('sweep_geometries15-Feb-2024-22-18-45.mat')
  
flux_data = [2.8, 1.9, 1.44, 1.039, 1.0, NaN];
flux_data_errors = [.175, .114, .067, .08466, .07237,NaN];

exportin = [1.246, 1.112, 1.013,1.057, 1.0, 1.0];

u0vals = [1.05, .925, .85, .8, .72, .75];

kmeasure = [2.11, 1.747, 1.412, 1.56, 1, 1];

kflat = 1.138;
Eflat = .7326;

Dfast = 1;
hmax=.1;
decay = 0;
box_height = 5;
box_width=15;

loss_full_simple = @(logTheta) fit_full_simple(exp(logTheta(1)),exp(logTheta(2)),Dfast, decay, u0vals, exportin, Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure,  kflat,  flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

loss_full_simple_nok = @(logTheta) fit_full_simple_nok(exp(logTheta(1)),exp(logTheta(2)),Dfast, decay, u0vals, exportin,  Eflat,  pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure,  kflat,  flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

loss_full = @(logTheta) fit_full(exp(logTheta(1)),exp(logTheta(2)),Dfast, decay, u0vals, exportin,  Eflat,pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure,  kflat,  flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

loss_nok = @(logTheta) fit_nok(exp(logTheta(1)),exp(logTheta(2)),Dfast, decay, u0vals, exportin,  Eflat,pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure,  kflat,  flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

loss_noexp = @(logTheta) fit_noexp(exp(logTheta(1)),exp(logTheta(2)),Dfast, decay, u0vals, exportin,  Eflat,pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure,  kflat,  flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

loss_noexpk = @(logTheta) fit_noexpk(exp(logTheta(1)),exp(logTheta(2)),Dfast, decay, u0vals, exportin,   Eflat,pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );


lb = [-4, -4];
ub = [3+log(Dfast), log(Dfast)];
init = [log(Dfast/10),log(Dfast/10)];

options = optimset('Display','iter','TolX',1e-3,'TolFun',1e-3);

[x_full_simple,fval_full_simple,exitflag_full_simple,output_full_simple] = fminsearchbnd(loss_full_simple,init,lb,ub,options);
[x_full_simple_nok,fval_full_simple_nok,exitflag_full_simple_nok,output_full_simple_nok] = fminsearchbnd(loss_full_simple_nok,init,lb,ub,options);

[x_full,fval_full,exitflag_full,output_full] = fminsearchbnd(loss_full,init,lb,ub,options);
[x_nok,fval_nok,exitflag_nok,output_nok] = fminsearchbnd(loss_nok,init,lb,ub,options);

[x_noexp,fval_noexp,exitflag_noexp,output_noexp] = fminsearchbnd(loss_noexp,init,lb,ub,options);
[x_noexpk,fval_noexpk,exitflag_noexpk,output_noexpk] = fminsearchbnd(loss_noexpk,init,lb,ub,options);

save('PDE_fit_save');


[~,out_full_simple] = fit_full_simple(exp(x_full_simple(1)),exp(x_full_simple(2)),Dfast, decay, u0vals, exportin,Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure,kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

[~,out_full_simple_nok] = fit_full_simple_nok(exp(x_full_simple_nok(1)),exp(x_full_simple_nok(2)),Dfast, decay, u0vals, exportin, Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

[~,out_full] = fit_full(exp(x_full(1)),exp(x_full(2)),Dfast, decay, u0vals, exportin,  Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure,kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

[~,out_nok] = fit_nok(exp(x_nok(1)),exp(x_nok(2)),Dfast, decay, u0vals, exportin,  Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

[~,out_noexp] = fit_noexp(exp(x_noexp(1)),exp(x_noexp(2)),Dfast, decay, u0vals, exportin,  Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

[~,out_noexpk] = fit_noexpk(exp(x_noexpk(1)),exp(x_noexpk(2)),Dfast, decay, u0vals, exportin, Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax );

save('PDE_fit_save');



%%%
Dfast_scale = 10;


figure('position',[200, 200, 800, 800]);
subplot(2,2,1)
hold on;

plot(out_full(1:5),'x--','DisplayName','full','Color',[254,97,0]/256,'LineWidth',1.25,'MarkerSize',10);
plot(out_full_simple(1:5),'+-','DisplayName','approx','Color',[254,97,0]/256,'LineWidth',1.25,'MarkerSize',10);
plot(out_full_simple_nok(1:5),'o-','DisplayName','approx, NPC sat','Color',[100,143,255]/256,'LineWidth',1.25,'MarkerSize',10);
plot(out_nok(1:5),'o--','DisplayName','full, NPC sat','Color',[100,143,255]/256, 'LineWidth',1.25,'MarkerSize',10);
plot(out_noexp(1:5),'-p','DisplayName','full, exp sat','Color',[120,94,240]/256, 'LineWidth',1.25,'MarkerSize',10)
plot(out_noexpk(1:5),'-s','DisplayName','full, exp & NPC sat','Color',[220,38,127]/256, 'LineWidth',1.25,'MarkerSize',10)

legend;

errorbar(1:6,flux_data, flux_data_errors,'o','MarkerSize', 8,'MarkerEdgeColor',[25,25,25]/256,'DisplayName','RPL13 data','LineWidth',2,...
    'Color','black');

 xticks([1,2,3,4,5,6,7])
 labels = {"R150","R200","R300","R300","R350","R400"};
 xticklabels(labels)
 ylabel('normalized flux out')
 xlim([0.5 5.5])
 legend box off;
 set(gca,'FontSize',14)
set(gca,'LineWidth',1.25)
ylim([0.8 4]);


subplot(2,2,2)

plot(1, fval_full,'x--','DisplayName','full','Color',[254,97,0]/256,'LineWidth',1.5,'MarkerSize',12);
hold on;
plot(2, fval_full_simple,'+-','DisplayName','approx','Color',[254,97,0]/256,'LineWidth',1.5,'MarkerSize',12);
plot(3, fval_full_simple_nok,'o-','DisplayName','approx, NPC sat','Color',[100,143,255]/256,'LineWidth',1.5,'MarkerSize',12);
plot(4, fval_nok,'o--','DisplayName','full, NPC sat','Color',[100,143,255]/256, 'LineWidth',1.5,'MarkerSize',12);
plot(5, fval_noexp,'-p','DisplayName','full, exp sat','Color',[120,94,240]/256, 'LineWidth',1.5,'MarkerSize',12)
plot(6, fval_noexpk,'-s','DisplayName','full, exp & NPC sat','Color',[220,38,127]/256, 'LineWidth',1.5,'MarkerSize',12)

xticks([1 2 3 4 5 6])
labels={"full","approx","approx, NPC sat", "full, NPC sat", "full, exp sat", "full exp & NPC sat"};
xticklabels(labels)
xlim([0.5 6.5]);
% legend box off;
 set(gca,'FontSize',12)
set(gca,'LineWidth',1.25)
set(gca,'YScale','log')
ylabel('mean squared error (MSE)')
box off;
grid on;
ylim([1e-3 1e0])
%ylim([1e-2,1e2])





subplot(2,2,3)

plot(1, Dfast_scale*exp(x_full(1)),'x--','DisplayName','full','Color',[254,97,0]/256,'LineWidth',1.5,'MarkerSize',12);
hold on;
plot(2, Dfast_scale*exp(x_full_simple(1)),'+-','DisplayName','approx','Color',[254,97,0]/256,'LineWidth',1.5,'MarkerSize',12);
plot(3, Dfast_scale*exp(x_full_simple_nok(1)),'o-','DisplayName','approx, NPC sat','Color',[100,143,255]/256,'LineWidth',1.5,'MarkerSize',12);
plot(4, Dfast_scale*exp(x_nok(1)),'o--','DisplayName','full, NPC sat','Color',[100,143,255]/256, 'LineWidth',1.5,'MarkerSize',12);
plot(5, Dfast_scale*exp(x_noexp(1)),'-p','DisplayName','full, exp sat','Color',[120,94,240]/256, 'LineWidth',1.5,'MarkerSize',12)
plot(6, Dfast_scale*exp(x_noexpk(1)),'-s','DisplayName','full, exp & NPC sat','Color',[220,38,127]/256, 'LineWidth',1.5,'MarkerSize',12)

xticks([1 2 3 4 5 6])
labels={"full","approx","approx, NPC sat", "full, NPC sat", "full, exp sat", "full exp & NPC sat"};
xticklabels(labels)
xlim([0.5 6.5]);
% legend box off;
 set(gca,'FontSize',12)
set(gca,'LineWidth',1.25)
set(gca,'YScale','log')
ylabel('fit k [um/s]')
box off;
grid on;
ylim([1e-1*Dfast_scale,1e2*Dfast_scale])


subplot(2,2,4)

plot(1, Dfast_scale*exp(x_full(2)),'x--','DisplayName','full','Color',[254,97,0]/256,'LineWidth',1.5,'MarkerSize',12);
hold on;
plot(2, Dfast_scale*exp(x_full_simple(2)),'+-','DisplayName','approx','Color',[254,97,0]/256,'LineWidth',1.5,'MarkerSize',12);
plot(3, Dfast_scale*exp(x_full_simple_nok(2)),'o-','DisplayName','approx, NPC sat','Color',[100,143,255]/256,'LineWidth',1.5,'MarkerSize',12);
plot(4, Dfast_scale*exp(x_nok(2)),'o--','DisplayName','full, NPC sat','Color',[100,143,255]/256, 'LineWidth',1.5,'MarkerSize',12);
plot(5, Dfast_scale*exp(x_noexp(2)),'-p','DisplayName','full, exp sat','Color',[120,94,240]/256, 'LineWidth',1.5,'MarkerSize',12)
plot(6, Dfast_scale*exp(x_noexpk(2)),'-s','DisplayName','full, exp & NPC sat','Color',[220,38,127]/256, 'LineWidth',1.5,'MarkerSize',12)

xticks([1 2 3 4 5 6])
labels={"full","approx","approx, NPC sat", "full, NPC sat", "full, exp sat", "full exp & NPC sat"};
xticklabels(labels)
xlim([0.5 6.5]);
% legend box off;
 set(gca,'FontSize',12)
set(gca,'LineWidth',1.25)
set(gca,'YScale','log')
ylabel('fit D_{slow} [um/s]')
box off;
grid on;
ylim([1e-2*Dfast_scale,1e0*Dfast_scale])




function [loss,out_ij] = fit_full(k,Dslow,Dfast, decay, u0vals, exportin,Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom_change, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax)


flux_out = zeros(1,6);
parfor i = 1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    x_nuc = squeeze(xvals_chrom_change(i,:,:));
    k_i = k*kmeasure(i);
    k_flat = k*kflat;
    exportin_i = exportin(i);
    exp_flat = Eflat;

    flux_out_i = solve_pde_flux_out(box_width,box_height,pillar_rad,pillar_dist,...
        pillar_height,...
        x_nuc,u0,k_i,k_flat,decay,Dfast,Dslow, cthick,cthick_nuc, exportin_i, exp_flat, hmax,0);

    flux_out(i) = flux_out_i;
end


out_ij = flux_out;
out_ij = out_ij/out_ij(5);
loss = nanmean((out_ij(1:5) - flux_data(1:5)).^2);

end


function [loss,out_ij] = fit_full_simple(k,Dslow,Dfast, decay, u0vals, exportin,Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom_change, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax)

simple_model_flux = @(D, k, L, u0) D*k*u0/(D+k*L);

flux_out = zeros(1,6);
parfor i = 1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    x_nuc = squeeze(xvals_chrom_change(i,:,:));
    k_i = k*kmeasure(i);

    flux_out_i = simple_model_flux(Dslow, k_i, cthick, u0)...
        *exportin(i);
 

    flux_out(i) = flux_out_i;
end

out_ij = flux_out;
out_ij = out_ij/out_ij(5);
loss = nanmean((out_ij(1:5) - flux_data(1:5)).^2);
end


function [loss,out_ij] = fit_full_simple_nok(k,Dslow,Dfast, decay, u0vals, exportin,Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom_change, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax)

simple_model_flux = @(D, k, L, u0) D*k*u0/(D+k*L);

flux_out = zeros(1,6);
parfor i = 1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    x_nuc = squeeze(xvals_chrom_change(i,:,:));
    k_i = k;
    

    flux_out_i = simple_model_flux(Dslow, k_i, cthick, u0)...
        *exportin(i);

    flux_out(i) = flux_out_i;
end


out_ij = flux_out;
out_ij = out_ij/out_ij(5);
loss = nanmean((out_ij(1:5) - flux_data(1:5)).^2);
end

function [loss,out_ij] = fit_nok(k,Dslow,Dfast, decay, u0vals, exportin,Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom_change, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax)


flux_out = zeros(1,6);
parfor i = 1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    x_nuc = squeeze(xvals_chrom_change(i,:,:));
    k_i = k;
    k_flat = k;
    exportin_i = exportin(i);
    exp_flat = Eflat;

    flux_out_i = solve_pde_flux_out(box_width,box_height,pillar_rad,pillar_dist,...
        pillar_height,...
        x_nuc,u0,k_i,k_flat,decay,Dfast,Dslow, cthick,cthick_nuc, exportin_i, exp_flat, hmax,0);

    flux_out(i) = flux_out_i;
end


out_ij = flux_out;
out_ij = out_ij/out_ij(end-1);
loss = nanmean((out_ij(1:5) - flux_data(1:5)).^2);
end

function [loss,out_ij] = fit_noexp(k,Dslow,Dfast, decay, u0vals, exportin,Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom_change, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax)


flux_out = zeros(1,7);
parfor i = 1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    x_nuc = squeeze(xvals_chrom_change(i,:,:));
    k_i = k*kmeasure(i);
    k_flat = k*kflat;
    exportin_i = 1;
    exp_flat = 1;

    flux_out_i = solve_pde_flux_out(box_width,box_height,pillar_rad,pillar_dist,...
        pillar_height,...
        x_nuc,u0,k_i,k_flat,decay,Dfast,Dslow, cthick,cthick_nuc, exportin_i, exp_flat, hmax,0);

    flux_out(i) = flux_out_i;
end


out_ij = flux_out;
out_ij = out_ij/out_ij(5);
loss = nanmean((out_ij(1:5) - flux_data(1:5)).^2);
end


function [loss,out_ij] = fit_noexpk(k,Dslow,Dfast, decay, u0vals, exportin,Eflat, pillar_radii, chromatin_thick_pillar, chromatin_thick_nuc,...
    xvals_chrom_change, kmeasure, kflat, flux_data, box_width, box_height,  pillar_height, pillar_dist, hmax)


flux_out = zeros(1,6);
parfor i = 1:6
    u0 = u0vals(i);
    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_nuc = chromatin_thick_nuc;
    x_nuc = squeeze(xvals_chrom_change(i,:,:));
    k_i = k;
    k_flat = k;
    exportin_i = 1;
    exp_flat = 1;

    flux_out_i = solve_pde_flux_out(box_width,box_height,pillar_rad,pillar_dist,...
        pillar_height,...
        x_nuc,u0,k_i,k_flat,decay,Dfast,Dslow, cthick,cthick_nuc, exportin_i, exp_flat, hmax,0);

    flux_out(i) = flux_out_i;
end


out_ij = flux_out;
out_ij = out_ij/out_ij(5);
loss = nanmean((out_ij(1:5) - flux_data(1:5)).^2);
end


