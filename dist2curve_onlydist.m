function distances = dist2curve_onlydist(curvexy,mapxy)

if size(mapxy,2)==1
    mapxy = mapxy';
end 
npts = size(mapxy,1);


distances = zeros(1,npts);
for i =1 :npts 

diff = vecnorm(curvexy-mapxy(i,:),2,2);
[min_d, ~] = min(diff);

%xy = curvexy(min_i,:);
distances(i)=min_d;

end 

end 