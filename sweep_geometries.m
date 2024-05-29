clear all;
close all;

box_width = 10;
pillar_dist = 3;
pillar_height = 1.5;

box_height=15;

M=151;
tmax=500;
dt=.025;

fster=1;
farea=5;
ftens=5;
fdown=.01;


init_pos = [0,4]';

pillar_radii = [.150 .200 .250 .300 .350 .400];

chromatin_thick_pillar = [ 0.3640	0.4973	0.5446	0.7003	0.7089 0.7089];

chromatin_thick_nuc = 0.7425;

no_area = [ 12.56	11.22	10.82	10.33	8.794	9.539	];
no_rads = sqrt(no_area/pi);

xvals_chrom = zeros(6, M,2);

tic;
parfor i = 1:6
    pillar_rad_i = pillar_radii(i);
    chrom_i = chromatin_thick_pillar(i);
    no_rad_i = no_rads(i);

    [x_bdry,y_bdry] = make_tot_bdry(box_width,box_height, pillar_rad_i,pillar_dist,pillar_height);

    x_chrom = make_blob_geometry(x_bdry,y_bdry,M, no_rad_i,init_pos,tmax,dt, ...
        fdown, fster, chrom_i, chromatin_thick_nuc, farea,ftens);

 %   x_chrom_fix = make_blob_geometry(x_bdry,y_bdry,M, no_rad_i,init_pos,tmax,dt, ...
  %      fdown, fster, chrom_fix,farea,ftens);

    xvals_chrom(i,:,:) = x_chrom;
  %  xvals_chrom_fix(i,:,:) = x_chrom_fix;

end
toc;
nowtime = datestr(datetime);
nowtime=strrep(nowtime,' ','-');
nowtime=strrep(nowtime,':','-');


filename = strcat('sweep_geometries',nowtime,'.mat');
save(filename);

figure;

for i = 1:6
    x_chrom = squeeze(xvals_chrom(i,:,:));
   % x_chrom_fix = squeeze(xvals_chrom_fix(i,:,:));

    pillar_rad_i = pillar_radii(i);
    [x_bdry,y_bdry] = make_tot_bdry(box_width,box_height, pillar_rad_i,pillar_dist,pillar_height);


    subplot(1,7,i);
    hold on;
    plot(x_bdry,y_bdry,'LineWidth',2,'color',[231,41,138]/256)

    plot(x_chrom(:,1),x_chrom(:,2),'LineWidth',2,'color',[102,166,30]/256);

    set(gca,'LineWidth',1.5)
    set(gcf,'color','w');
    axis equal;
        ylim([0 6])
    xlim([-3, 3]+init_pos(1));
    % 
    % subplot(2,7,7+i);
    % hold on;
    % plot(x_bdry,y_bdry,'LineWidth',2,'color',[231,41,138]/256)
    % 
    % plot(x_chrom_fix(:,1),x_chrom_fix(:,2),'LineWidth',2,'color',[102,166,30]/256);
    % 
    % set(gca,'LineWidth',1.5)
    % set(gcf,'color','w');
    % axis equal;
    %     ylim([0 6])
    % xlim([-3, 3]);

end

%
% tic;
% x= make_blob_geometry(x_bdry,y_bdry,M, init_rad,init_pos,tmax,dt, ...
%     fdown, fster, lster,farea,ftens);
% toc;
% figure;
% hold on;
%             plot(x_bdry,y_bdry,'LineWidth',2,'color',[231,41,138]/256)
%
%             plot(x(:,1),x(:,2),'LineWidth',2,'color',[102,166,30]/256);
%             ylim([0 6])
%             xlim([-5, 5]);
%             set(gca,'LineWidth',1.5)
%             set(gcf,'color','w');
%             axis equal;
