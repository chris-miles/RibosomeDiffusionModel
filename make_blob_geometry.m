function x= make_blob_geometry(x_bdry,y_bdry,M, init_rad,init_pos,tmax,dt, ...
    fdown, fster, lster, lster_base, farea,ftens)


%M=151;
r = init_rad;

x = r*[cos(linspace(pi/2,pi/2+ 2*pi,M));sin(linspace(pi/2,pi/2+2*pi,M))]+init_pos;

x= x';

nT=round(tmax/dt);


%fdown = @(t).05;%*exp(-t/tmax);
%k=25;

%fster=1;
%lster = .25;
%l_aster_big = 2*lster;
%fster_func = @(d) fster*(d-lster)*(d<=l_aster_big);
%lster = .15;

fster_func = @(d) fster*(1-d/lster)*(d<=lster);
fster_func_base = @(d) fster*(1-d/lster_base)*(d<=lster_base);


[xarc,~] = arclength(x(:,1),x(:,2));
target_dr = xarc/M;%2*pi*r/(M-1);
target_area = pi*(r^2);

%farea=5;


t = 0;
x = flipud(x);


for i = 1:nT
    t = t+dt;

    %disp(num2str(t/tmax));


    x(M,:) = x(1,:);
    [xarc,sarc] = arclength(x(:,1),x(:,2));
    xnew = interp1([0;cumsum(sarc)],x,linspace(0,xarc,M));
    x=xnew;

    x(M,:) = x(1,:);
    x = flipud(x);


    diff_x = diff(x);%gradient(x')';
    tensmag = ftens*(vecnorm(diff_x,2,2)-target_dr).*(diff_x).*(vecnorm(diff_x,2,2).^(-1));
    Ftens = zeros(M,2);
    for m = 2:M-1
        Ftens(m,:) = tensmag(m-1,:) - tensmag(m,:);
    end

    Ftens(1,:) = -tensmag(1,:) + tensmag(end,:);
    Ftens(M,:) = Ftens(1,:);
    Ftens=-Ftens;


    % enforce symmetry

    %Ftens(abs(Ftens)<1e-8) =0;

    [~,normalvec,~,~,~] = frenet(x(:,1),x(:,2));


    Fdown = -[0,fdown];

    %F_area = farea*((polyarea(x(:,1),x(:,2)) - target_area))*normalvec(:,1:2);%farea*log(polyarea(x(:,1),x(:,2))/target_area)*normalvec(:,1:2);
    F_area = farea*log(polyarea(x(:,1),x(:,2))/target_area)*normalvec(:,1:2);


    F_ster = zeros(M,2);
    for m = 1:M
        x_m = x(m,:);
        [closest_pt,dist_to_surf] = dist2curve([x_bdry;y_bdry]',x_m);

        force_dir = -(closest_pt-x_m)/vecnorm(closest_pt-x_m);

        force_mag = fster_func(dist_to_surf);
        force_mag_base = fster_func_base(x_m(2));
        
        if force_mag_base > force_mag
            force_mag = force_mag_base;
            force_dir = [0,1];
        end 

        F_ster(m,:) = force_mag*force_dir;
    end

    Ftot = Fdown+Ftens+F_ster+F_area;

    x = x+dt*Ftot;


end


end