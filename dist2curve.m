function [xy,distance] = dist2curve(curvexy,mapxy)

diff = vecnorm(curvexy-mapxy,2,2);
[min_d, min_i] = min(diff);

xy = curvexy(min_i,:);
distance=min_d;

end 