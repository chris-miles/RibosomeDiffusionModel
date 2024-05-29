function [x_bdry,y_bdry,labels] = make_tot_bdry(box_width,box_height, pillar_rad,pillar_dist,pillar_height)

h = pillar_height;
if pillar_rad>0
    dx = .01;


    lower_bdry1_x = -box_width/2:dx:(-pillar_rad-pillar_dist);
    lower_bdry1_y = zeros(size(lower_bdry1_x));

    lower_bdry2_y = 0:dx:h;
    lower_bdry2_x = (-pillar_rad-pillar_dist)*ones(size(lower_bdry2_y));

    lower_bdry3_x = (-pillar_rad-pillar_dist):dx:(+pillar_rad-pillar_dist);
    lower_bdry3_y = h*ones(size(lower_bdry3_x));

    lower_bdry4_y = flip(0:dx:h);
    lower_bdry4_x = (+pillar_rad-pillar_dist)*ones(size(lower_bdry4_y));

    lower_bdry5_x = (+pillar_rad-pillar_dist):dx:(-pillar_rad);
    lower_bdry5_y = zeros(size(lower_bdry5_x));

    lower_bdry6_y = 0:dx:h;
    lower_bdry6_x = (-pillar_rad)*ones(size(lower_bdry6_y));

    lower_bdry7_x = (-pillar_rad):dx:(+pillar_rad);
    lower_bdry7_y = h*ones(size(lower_bdry7_x));

    lower_bdry8_y = flip(0:dx:h);
    lower_bdry8_x = (+pillar_rad)*ones(size(lower_bdry8_y));

    lower_bdry9_x = (pillar_rad):dx:(-pillar_rad+pillar_dist);
    lower_bdry9_y = zeros(size(lower_bdry9_x));

    lower_bdry10_y = 0:dx:h;
    lower_bdry10_x = (-pillar_rad+pillar_dist)*ones(size(lower_bdry10_y));

    lower_bdry11_x =  (-pillar_rad+pillar_dist):dx: (+pillar_rad+pillar_dist);
    lower_bdry11_y = h*ones(size(lower_bdry11_x));

    lower_bdry12_y = flip(0:dx:h);
    lower_bdry12_x = (+pillar_rad+pillar_dist)*ones(size(lower_bdry12_y));

    lower_bdry13_x = (+pillar_rad+pillar_dist):dx:box_width/2;
    lower_bdry13_y = zeros(size(lower_bdry13_x));

    lower_bdry14_y = 0:dx:box_height;
    lower_bdry14_x = (box_width/2)*ones(size(lower_bdry14_y));

    lower_bdry15_x = flip(-box_width/2:dx:box_width/2);
    lower_bdry15_y = box_height*ones(size(lower_bdry15_x));

    lower_bdry16_y = flip(0:dx:box_height);
    lower_bdry16_x = (-box_width/2)*ones(size(lower_bdry16_y));


    x_bdry = [lower_bdry1_x,lower_bdry2_x,lower_bdry3_x,lower_bdry4_x,lower_bdry5_x,lower_bdry6_x, lower_bdry7_x, lower_bdry8_x, ...
        lower_bdry9_x, lower_bdry10_x,lower_bdry11_x,lower_bdry12_x,lower_bdry13_x,lower_bdry14_x,lower_bdry15_x,lower_bdry16_x];

    y_bdry = [lower_bdry1_y,lower_bdry2_y,lower_bdry3_y,lower_bdry4_y,lower_bdry5_y,lower_bdry6_y, lower_bdry7_y, lower_bdry8_y, ...
        lower_bdry9_y, lower_bdry10_y,lower_bdry11_y,lower_bdry12_y,lower_bdry13_y,lower_bdry14_y,lower_bdry15_y,lower_bdry16_y];


    labels = [1*ones(size(lower_bdry1_x)), 2*ones(size(lower_bdry2_x)),3*ones(size(lower_bdry3_x)),4*ones(size(lower_bdry4_x)),5*ones(size(lower_bdry5_x)),...
      6*ones(size(lower_bdry6_x)), 7*ones(size(lower_bdry7_x)),8*ones(size(lower_bdry8_x)), 9*ones(size(lower_bdry9_x)), 10*ones(size(lower_bdry10_x)),...
      11*ones(size(lower_bdry11_x)), 12*ones(size(lower_bdry12_x)),13*ones(size(lower_bdry13_x)),14*ones(size(lower_bdry14_x)),15*ones(size(lower_bdry15_x)),...
      16*ones(size(lower_bdry16_x))]';

else
    dx = .01;%box_width/500;

    x1_bdry = -box_width/2:dx:box_width/2;
    y1_bdry = zeros(size(x1_bdry));

    lower_bdry14_y = 0:dx:box_height;
    lower_bdry14_x = (box_width/2)*ones(size(lower_bdry14_y));

    lower_bdry15_x = flip(-box_width/2:dx:box_width/2);
    lower_bdry15_y = (box_height)*ones(size(lower_bdry15_x));

    lower_bdry16_y = flip(0:dx:box_height);
    lower_bdry16_x = (-box_width/2)*ones(size(lower_bdry16_y));


    x_bdry = [x1_bdry, ...
        lower_bdry14_x ,lower_bdry15_x ,lower_bdry16_x ];

    y_bdry = [y1_bdry,lower_bdry14_y ,lower_bdry15_y ,lower_bdry16_y ];
    labels = [ones(size(x1_bdry)),2*ones(size(lower_bdry14_x)),3*ones(size(lower_bdry15_x)),...
        4*ones(size(lower_bdry16_x))]';

end
end
