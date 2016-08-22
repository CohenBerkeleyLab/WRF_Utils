function [ xloncorn, xlatcorn ] = wrf_grid_corners( xlon, xlat )
%WRF_GRID_CORNERS Computes corners of WRF grid cells.
%   [ XLONCORN, XLATCORN ] = WRF_GRID_CORNERS( XLON, XLAT )
%       Takes as input the grid center points XLON and XLAT. Returns the
%       4-by-m-by-n arrays XLONCORN and XLATCORN where m and n are the
%       first and second dimension lengths of XLON and XLAT.


xloncorn = nan([4, size(xlon)]);
xlatcorn = nan([4, size(xlat)]);
% The interior corners will be shared amongst neighbor points, so we only
% need to calculate them once.
interior_loncorn = (xlon(1:end-1,1:end-1)+xlon(2:end,1:end-1)+xlon(2:end,2:end)+xlon(1:end-1,2:end))/4;
interior_latcorn = (xlat(1:end-1,1:end-1)+xlat(2:end,1:end-1)+xlat(2:end,2:end)+xlat(1:end-1,2:end))/4;

xloncorn(1,2:end,2:end) = interior_loncorn;
xlatcorn(1,2:end,2:end) = interior_latcorn;

xloncorn(2,1:end-1,2:end) = interior_loncorn;
xlatcorn(2,1:end-1,2:end) = interior_latcorn;

xloncorn(3,1:end-1,1:end-1) = interior_loncorn;
xlatcorn(3,1:end-1,1:end-1) = interior_latcorn;

xloncorn(4,2:end,1:end-1) = interior_loncorn;
xlatcorn(4,2:end,1:end-1) = interior_latcorn;

% Now we do the edges, ignoring the extreme corners b/c we don't actually 
% have enough information right now. Assume that the spacing between the
% first and second column is the same as between the second and third, and
% so on.
xloncorn(2,1:end-1,1) = compute_by_difference( xloncorn(2,1:end-1,2), xloncorn(2,1:end-1,3) );
xloncorn(1,2:end,1) = xloncorn(2,1:end-1,1);
xlatcorn(2,1:end-1,1) = compute_by_difference( xlatcorn(2,1:end-1,2), xlatcorn(2,1:end-1,3) );
xlatcorn(1,2:end,1) = xlatcorn(2,1:end-1,1);

xloncorn(3,1:end-1,end) = compute_by_difference( xloncorn(3,1:end-1,end-1), xloncorn(3,1:end-1,end-2) );
xloncorn(4,2:end,end) = xloncorn(3,1:end-1,end);
xlatcorn(3,1:end-1,end) = compute_by_difference( xlatcorn(3,1:end-1,end-1), xlatcorn(3,1:end-1,end-2) );
xlatcorn(4,2:end,end) = xlatcorn(3,1:end-1,end);

xloncorn(4,1,1:end-1) = compute_by_difference( xloncorn(4,2,1:end-1), xloncorn(4,3,1:end-1) );
xloncorn(1,1,2:end) = xloncorn(4,1,1:end-1);
xlatcorn(4,1,1:end-1) = compute_by_difference( xlatcorn(4,2,1:end-1), xlatcorn(4,3,1:end-1) );
xlatcorn(1,1,2:end) = xlatcorn(4,1,1:end-1);

xloncorn(3,end,1:end-1) = compute_by_difference( xloncorn(3,end-1,1:end-1), xloncorn(3,end-2,1:end-1) );
xloncorn(2,end,2:end) = xloncorn(3,end,1:end-1);
xlatcorn(3,end,1:end-1) = compute_by_difference( xlatcorn(3,end-2,1:end-1), xlatcorn(3,end-2,1:end-1) );
xlatcorn(2,end,2:end) = xlatcorn(3,end,1:end-1);

% Finally the four remaining corners, now that we've filled in the middles
% of the edges, we have what we need. We'll use the average from both
% directions to determine this.
xloncorn(1,1,1) = mean([compute_by_difference(xloncorn(1,2,1), xloncorn(1,3,1)), compute_by_difference(xloncorn(1,1,2), xloncorn(1,1,3))]);
xlatcorn(1,1,1) = mean([compute_by_difference(xlatcorn(1,2,1), xlatcorn(1,3,1)), compute_by_difference(xlatcorn(1,1,2), xlatcorn(1,1,3))]);

xloncorn(2,end,1) = mean([compute_by_difference(xloncorn(2,end-1,1), xloncorn(2,end-2,1)), compute_by_difference(xloncorn(2,end,2), xloncorn(2,end,3))]);
xlatcorn(2,end,1) = mean([compute_by_difference(xlatcorn(2,end-1,1), xlatcorn(2,end-2,1)), compute_by_difference(xlatcorn(2,end,2), xlatcorn(2,end,3))]);

xloncorn(3,end,end) = mean([compute_by_difference(xloncorn(3,end-1,end), xloncorn(3,end-2,end)), compute_by_difference(xloncorn(3,end,end-1), xloncorn(3,end,end-2))]);
xlatcorn(3,end,end) = mean([compute_by_difference(xlatcorn(3,end-1,end), xlatcorn(3,end-2,end)), compute_by_difference(xlatcorn(3,end,end-1), xlatcorn(3,end,end-2))]);

xloncorn(4,1,end) = mean([compute_by_difference(xloncorn(4,2,end), xloncorn(4,3,end)), compute_by_difference(xloncorn(4,1,end-1), xloncorn(4,1,end-2))]);
xlatcorn(4,1,end) = mean([compute_by_difference(xlatcorn(4,2,end), xlatcorn(4,3,end)), compute_by_difference(xlatcorn(4,1,end-1), xlatcorn(4,1,end-2))]);
end

function [ out ] = compute_by_difference(nearest, next_nearest)
out = nearest - (next_nearest - nearest);
end