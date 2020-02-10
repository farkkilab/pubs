function density = densityfield(xydata,pospop)

k = 10; %this parameter determines how many neighbors to use to estimate density

%density is estimated as proportional to 1/distance^2, i.e. if you're
%dense, your k'th neighbor is still VERY close. It can be proved that the
%expectation value of density ought to depend on the k'th neighbors
%distance in this fashion

tic
[~,distall] = knnsearch(xydata,xydata,'k',k,'NSMethod','kdtree');
basaldensity = 1./(distall(:,k)).^2; %calculate density of all cells, to normalize later (loss of cells, or smaller cells become packed more densely)
[~,distpop] = knnsearch(xydata(pospop,:),xydata,'k',k,'NSMethod','kdtree');
popdensity = 1./(distpop(:,k)).^2; %calculate density of the chosen population, for each spatial position (i.e. every cell's position)
density = popdensity./basaldensity; %normalize
toc
end