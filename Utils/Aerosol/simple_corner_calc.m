function [ modis_latcorn, modis_loncorn ] = simple_corner_calc( modis_lat, modis_lon )
%Extremely simple calculation of pixel corners; not recommended for
%publishing but probably okay for pre-gridded data

if ~all(size(modis_lat)==size(modis_lon)); error('simple_corner_calc:mismatch','Input lat and lon matrices are not the same size'); end

modis_latcorn = zeros(4,size(modis_lat,1),size(modis_lat,2));
modis_loncorn = zeros(4,size(modis_lat,1),size(modis_lat,2));

%1 will be the upper left corner, 2 the upper right, 3 the lower right, 4
%the lower left

modis_latcorn(1,:,:) = [modis_lat(1,:) + (modis_lat(1,:)-modis_lat(2,:))/2 ; modis_lat(2:end,:) + (modis_lat(1:(end-1),:) - modis_lat(2:end,:))./2];
modis_latcorn(2,:,:) = [modis_lat(1,:) + (modis_lat(1,:)-modis_lat(2,:))/2 ; modis_lat(2:end,:) + (modis_lat(1:(end-1),:) - modis_lat(2:end,:))./2];
modis_latcorn(3,:,:) = [modis_lat(1:(end-1),:) - (modis_lat(1:(end-1),:) - modis_lat(2:end,:))./2; modis_lat(end,:) - (modis_lat((end-1),:)-modis_lat(end,:))/2 ];
modis_latcorn(4,:,:) = [modis_lat(1:(end-1),:) - (modis_lat(1:(end-1),:) - modis_lat(2:end,:))./2; modis_lat(end,:) - (modis_lat((end-1),:)-modis_lat(end,:))/2 ];

modis_loncorn(1,:,:) = [modis_lon(:,1) + (modis_lon(:,1)-modis_lon(:,2))/2, modis_lon(:,2:end) + (modis_lon(:,1:(end-1)) - modis_lon(:,2:end))./2];
modis_loncorn(2,:,:) = [modis_lon(:,1:(end-1)) - (modis_lon(:,1:(end-1))-modis_lon(:,2:end))/2, modis_lon(:,end) - (modis_lon(:,(end-1)) - modis_lon(:,end))./2];
modis_loncorn(3,:,:) = [modis_lon(:,1:(end-1)) - (modis_lon(:,1:(end-1))-modis_lon(:,2:end))/2, modis_lon(:,end) - (modis_lon(:,(end-1)) - modis_lon(:,end))./2];
modis_loncorn(4,:,:) = [modis_lon(:,1) + (modis_lon(:,1)-modis_lon(:,2))/2, modis_lon(:,2:end) + (modis_lon(:,1:(end-1)) - modis_lon(:,2:end))./2];

end

