clear all;
%close all;

load('sweep_geometries15-Feb-2024-22-18-45.mat')
 set(gcf, 'defaultFigureRenderer', 'painters')

figure('position',[200,0, 400,1200]);

set(gcf, 'defaultFigureRenderer', 'painters')
labels = {"R150","R200","R300","R300","R350","R400"};
box_width = 15;
box_height=5;

for i = 1:6


    pillar_rad = pillar_radii(i);
    cthick = chromatin_thick_pillar(i);
    cthick_bottom = chromatin_thick_nuc;

    [xb_tot, yb_tot,b_labels] = make_tot_bdry(box_width,box_height,pillar_rad,pillar_dist,pillar_height);

    x_nuc = squeeze(xvals_chrom(i,:,:));

    if pillar_rad > 0

        rect_out = [3 % indicates rectangle
            4 % line segments, always 4
            -box_width/2 % x0
            box_width/2  % x1
            box_width/2  % x2
            -box_width/2 % x3
            0 % y0
            0 % y1
            box_height % y2
            box_height]; %y3


        pillar1 = [3 % indicates rectangle
            4 % line segments, always 4
            -pillar_dist-pillar_rad % x0
            -pillar_dist-pillar_rad  % x1
            -pillar_dist+pillar_rad  % x2
            -pillar_dist+pillar_rad % x3
            0 % y0
            pillar_height % y1
            pillar_height % y2
            0]; %y3

        pillar2 = [3 % indicates rectangle
            4 % line segments, always 4
            -pillar_rad % x0
            -pillar_rad  % x1
            +pillar_rad  % x2
            +pillar_rad % x3
            0 % y0
            pillar_height % y1
            pillar_height % y2
            0]; %y3

        pillar3 = [3 % indicates rectangle
            4 % line segments, always 4
            pillar_dist-pillar_rad % x0
            pillar_dist-pillar_rad  % x1
            pillar_dist+pillar_rad  % x2
            pillar_dist+pillar_rad % x3
            0 % y0
            pillar_height % y1
            pillar_height % y2
            0]; %y3

        p1 = [2
            M-1
            x_nuc(1:M-1,1)
            x_nuc(1:M-1,2)];



        rect_out = [rect_out;zeros(length(p1) - length(rect_out),1)];
        pillar1 = [pillar1;zeros(length(p1) - length(pillar1),1)];
        pillar2 = [pillar2;zeros(length(p1) - length(pillar2),1)];
        pillar3 = [pillar3;zeros(length(p1) - length(pillar3),1)];


        gd = [p1,pillar1,pillar2,pillar3,rect_out];

        ns = char('p1','pillar1','pillar2','pillar3','rect_out');
        ns = ns';

        sf = 'rect_out-pillar1-pillar2-pillar3-p1';

        [dl,~] = decsg(gd,sf,ns);

    else


        rect_out = [3 % indicates rectangle
            4 % line segments, always 4
            -box_width/2 % x0
            box_width/2  % x1
            box_width/2  % x2
            -box_width/2 % x3
            0 % y0
            0 % y1
            box_height % y2
            box_height]; %y3

        p1 = [2
            M-1
            x_nuc(1:M-1,1)
            x_nuc(1:M-1,2)];


        rect_out = [rect_out;zeros(length(p1) - length(rect_out),1)];

        gd = [p1,rect_out];

        ns = char('p1','rect_out');
        ns = ns';

        sf = 'rect_out-p1';

        [dl,~] = decsg(gd,sf,ns);

    end

    %pdegplot(dl,"EdgeLabels","on","FaceLabels","on")
    %axis equal



    model = createpde;
    geometryFromEdges(model,dl);

    nEdges = model.Geometry.NumEdges;

    reflect_edges = [M-1:nEdges-1];
    nuc_edges = [1:M-1,nEdges];

    [x_bdry,y_bdry] = make_tot_bdry(box_width, box_height,pillar_rad,pillar_dist,pillar_height);


    dist_to_bdry = @(x) dist2curve_onlydist([x_bdry;y_bdry]',x);
    dist_to_bottom = @(x) x(2);
    dist_to_nuc = @(x) dist2curve_onlydist([x_nuc(:,1),x_nuc(:,2)],x);
   


    %dist_to_bdry = @(x) dist2curve_onlydist([xb_tot;yb_tot]',x);


    msh=generateMesh(model,"Hmax",.25);

    nMesh = length(msh.Nodes);
    distBank = zeros(nMesh,1);
    distBank2 = zeros(nMesh,1);
    distBank3 = zeros(nMesh,1);

    for j = 1:nMesh
        distBank(j) = dist_to_bdry(msh.Nodes(:,j));
        distBank2(j) = dist_to_bottom(msh.Nodes(:,j));
         distBank3(j) = dist_to_nuc(msh.Nodes(:,j));       
    end

    hetero_region = (distBank<cthick_bottom)|(distBank2<cthick_bottom)|(distBank3<cthick);

    pillar_rad_i = pillar_radii(i);


    subplot(6,1,i);
    hold on;
    scatter(msh.Nodes(1,hetero_region)-init_pos(1),msh.Nodes(2,hetero_region),5,'filled','MarkerFaceColor',[102,166,30]/256,...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
    plot(x_bdry-init_pos(1),y_bdry,'LineWidth',2,'color',[55,126,184]/256 )

    ppp=plot(polyshape(x_nuc(:,1)-init_pos(1),x_nuc(:,2)));%,'LineWidth',2,'color',[231,41,138]/256);
    ppp.EdgeAlpha=0;
    ppp.FaceColor=[231,41,138]/256;
    ppp.FaceAlpha=1;

    %for j = 1:nMesh
    % if distBank(j) < cthick


    %   else
    %        scatter(msh.Nodes(1,j),msh.Nodes(2,j),5,'filled','MarkerFaceColor',[202,166,30]/256,...
    %  'MarkerEdgeColor','none');
    % end

    % end
    set(gca,'LineWidth',1.5)
    set(gcf,'color','w');
    axis equal;
    ylim([0 6])
    xlim([-8, 8]);
    title(labels{i},'FontWeight','normal');
set(gcf, 'defaultFigureRenderer', 'painters')
axis off;

end







