function [tot_sol] = solve_pde_flux_out(box_width,box_height, ...
    pillar_rad,pillar_dist,pillar_height,...
    x_nuc,u0,kpillar,kflat, decay,Dfast,Dslow, cthick,cthick_surf,epillar, eflat,hmax,plotflag)

%set(gcf, 'defaultFigureRenderer', 'painters')

if nargin == 17
    plotflag=0;
end 


[xb_tot, yb_tot,b_labels] = make_tot_bdry(box_width,box_height,pillar_rad,...
    pillar_dist,pillar_height);

[xb_box_only, yb_box_only] = make_tot_bdry(box_width,box_height,0,...
    pillar_dist,pillar_height);


M = length(x_nuc);

if pillar_rad > 0


vert_b = (mod(b_labels,2)==0)&b_labels<14;
horz_b = (mod(b_labels,2)==1)&b_labels<14;

vert_t = (mod(b_labels,2)==0)&b_labels>=14;
horz_t = (mod(b_labels,2)==1)&b_labels>=14;

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


vert_b = (mod(b_labels,2)==0)&b_labels<=1;
horz_b = (mod(b_labels,2)==1)&b_labels<=1;

vert_t = (mod(b_labels,2)==0)&b_labels>1;
horz_t = (mod(b_labels,2)==1)&b_labels>1;



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

% pillars are 150, 151, 152, 153, 154 SIDES
%              158 159 166 tops


%[150, 151, 152, 153, 154, 158, 159, 166];


model = createpde;
geometryFromEdges(model,dl);

nEdges = model.Geometry.NumEdges;


if pillar_rad>0
reflect_edges = [M-1:nEdges-2,nEdges];
nuc_edges = [1:M-2,nEdges-1];
else
reflect_edges = [M-1:nEdges-1];
nuc_edges = [1:M-2,nEdges];
end 

pillar_edges = [M-1:M+3, M+7, M+8, nEdges-2,nEdges];

reflect_minus_pillar = setdiff(reflect_edges, pillar_edges);

%u0=1;
%k=100;
%decay=.1;

%Dfast=100;
%Dslow=100;
ss=10;


applyBoundaryCondition(model,"dirichlet", ...
    "Edge",nuc_edges,"u",u0);

% j + k*u = 0

applyBoundaryCondition(model,"neumann", ...
    "Edge",reflect_minus_pillar,"g",0,"q",kflat);


applyBoundaryCondition(model,"neumann", ...
    "Edge",pillar_edges,"g",0,"q",kpillar);


%applyBoundaryCondition(model,"dirichlet", ...
%                             "Edge",reflect_edges,"u",0);


%[x_bdry,y_bdry] = make_bdry(box_width,pillar_rad,pillar_dist,pillar_height);

dist_to_bdry = @(x) dist2curve_onlydist([xb_tot;yb_tot]',x);
dist_to_box = @(x) dist2curve_onlydist([xb_box_only;yb_box_only]',x);
dist_to_nuc = @(x) dist2curve_onlydist([x_nuc(:,1),x_nuc(:,2)],x);


%dist_to_bdry = @(x) dist2curve_onlydist([xb_tot;yb_tot]',x);


msh=generateMesh(model,"Hmax",hmax);

nMesh = length(msh.Nodes);
distBank = zeros(nMesh,1);

for i = 1:nMesh
    distToPillar = dist_to_bdry(msh.Nodes(:,i))-cthick;
    distToBox = dist_to_box(msh.Nodes(:,i))-cthick_surf;
    distToNuc = dist_to_nuc(msh.Nodes(:,i)) -cthick;
    if distToPillar< 0
        distToPillar=0;
    end
    if distToBox<0
        distToBox=0;
    end 
    if distToNuc<0
        distToNuc=0;
    end 
    
    distBank(i) = min([distToPillar,distToBox,distToNuc]); 

end


interped_dist_to_bdry = scatteredInterpolant(msh.Nodes',distBank);%@(x) interp2(msh.Nodes(1,:),msh.Nodes(2,:),distbank,x(1),x(2));

Dfuncx = @(x) Dslow+(Dfast-Dslow)*(tanh(ss*(x))); 
%Dfuncx = @(x) Dslow+(Dfast-Dslow)*(tanh(ss*(x)))+(Dfast-Dslow)*exp(-ss*100*x);

Dfuncloc = @(location,state) Dfuncx(interped_dist_to_bdry([location.x;location.y]'))';


specifyCoefficients(model,"m",0,"d",0,"c",Dfuncloc,"a",decay,"f",0);

results = solvepde(model);
u = results.NodalSolution;

[cgradx,cgrady] = evaluateCGradient(results);


pillar_labels = [2,3,4,6,7,8,10,11,12];

pillar_boundary_indicies = ismember(b_labels,pillar_labels);
flat_bdry_indicies = setdiff(1:length(b_labels), find(pillar_boundary_indicies));

interped_sol = interpolateSolution(results, xb_tot,yb_tot);

u_on_pillar = interped_sol(pillar_boundary_indicies);
u_flat = interped_sol(flat_bdry_indicies);

dx_bdry = xb_tot(2)-xb_tot(1);


tot_sol = sum(u_on_pillar)*kpillar*epillar*dx_bdry + ...
    sum(u_flat)*kflat*eflat*dx_bdry;



%axis equal


%cmag = vecnorm([cgradx,cgrady],2,2);

if plotflag 

maxtot = 1.3799;
mintot=.4896;


addpath('cbrewer')
%colormap3 = flipud(cbrewer('seq','YlGnBu',100));

%colormap3=viridis(100);
colormap3=parula(100);
%colormap3=cmasher('emerald',100);
%colormap3 = flipud(cbrewer('seq','PuRd',150,'linear'));
%colormap3 = colormap3(25:125,:);%flipud(cbrewer('seq','PuRd',100,'linear'));

scaled_tot_sol = (tot_sol-mintot)/(maxtot-mintot);
color3 = abs(colormap3(round(99*scaled_tot_sol)+1,:));



figure('position',[200,200,1500,600],'renderer','painters'); 
hold on;
fill(.6*[-box_width,box_width,box_width,-box_width], ...
    [-0.2*box_height,-0.2*box_height,1.2*box_height,1.2*box_height],color3,'LineStyle','none');

%subplot(1,2,1)

title('concentration','FontWeight','normal')
uplot=pdeplot(model,"XYData",u,"Contour","on");
hold on;

%for jjj=2:(length(uplot)-1)
%uplot(jjj).Color='w';
%end

axis equal;
xlim([-box_width/2,box_width/2])
ylim([0,box_height])

colormap1=viridis(100);
%colormap1 = cmasher('bubblegum',100);
%colormap1=flipud(cbrewer('seq','YlGnBu',100,'linear'));
colormap(colormap1);

clim([0 1])
colorbar off
axis off;
colormap2=plasma(100);

%colormap2=turbo(100);
%colormap2 = cmasher('neon',100);
%colormap2=cbrewer('seq','RdPu',100,'linear');

max_val = log(.2944); % hard coded from r=150 sim max flux.
min_val = log(7e-5);



color_interp = @(x) round(99*(log(x)-min_val)/(max_val-min_val))+1;


pillar_colors = colormap2(color_interp(kpillar*epillar*interped_sol(pillar_boundary_indicies)),:);
flat_colors = colormap2(color_interp(kflat*eflat*interped_sol(flat_bdry_indicies)),:);


scatter3(xb_tot(pillar_boundary_indicies),yb_tot(pillar_boundary_indicies), ...
    kpillar*epillar*interped_sol(pillar_boundary_indicies),85, pillar_colors, ...
    'filled',"marker","square")
hold on
scatter3(xb_tot(flat_bdry_indicies),yb_tot(flat_bdry_indicies), ...
    kflat*eflat*interped_sol(flat_bdry_indicies),85, flat_colors, ...
  'filled',"marker","square")
view(2);


%xlim([-2*pillar_dist, 2*pillar_dist])
%ylim([0, 5])

% 
% subplot(1,2,2)
% title('flux','FontWeight','normal')
% 
% cplot=pdeplot(model,"XYData",cmag,"Contour","off","FlowData",-[cgradx,cgrady]);
% %cplot(2).Color=[252,141,98]/256;
% cplot(2).Color='w';
% 
% axis equal;
% xlim([-box_width/2,box_width/2])
% ylim([0,box_height])
% colormap plasma;
% axis off;
% %xlim([-2*pillar_dist, 2*pillar_dist])
% %ylim([0, 5])
% end 

%[cgrad_bdry_x, cgrad_bdry_y] = evaluateCGradient(results,xb_tot,yb_tot);
%dx_bdry = xb_tot(2)-xb_tot(1);

%vert_tot_b = sum(abs(cgrad_bdry_x(vert_b)))*dx_bdry;
%horz_tot_b = sum(abs(cgrad_bdry_y(horz_b)))*dx_bdry;

%vert_tot_t = sum(abs(cgrad_bdry_x(vert_t)))*dx_bdry;
%horz_tot_t = sum(abs(cgrad_bdry_y(horz_t)))*dx_bdry;

%flux_out = [vert_tot_b, horz_tot_b,vert_tot_t, horz_tot_t];


%set(gcf, 'defaultFigureRenderer', 'painters');
%axis off;

%sum_flux_out = sum(flux_out)
addpath('cbrewer')

v_colormap = colormap1;%cmasher('lavender',100);%viridis(100);

fill(x_nuc(:,1),x_nuc(:,2),v_colormap(round(u0*99/1.05)+1,:))


end 