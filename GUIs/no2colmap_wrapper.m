%no2colmap_wrapper
% This is a wrapper script for the no2_column_map_GUI function; it takes
% the output from the GUI script and parses it into the no2_column_map_2014
% function.  This allows the function to execute in the command window, so
% that the user can monitor its progress through printed messages and
% cancel the script with Control+C if necessary.  This will bring the output
% variables from the function into the workspace as well.
%
%   Josh Laughner <joshlaugh5@gmail.com>

output = no2_column_map_GUI;
[cb, NO2_GRID, LON_GRID, LAT_GRID] = no2_column_map_2014(output.start_date, output.end_date, output.lonbdy, output.latbdy,...
    'mapfield',output.data_field,'resolution',output.grid_resolution,'projection',output.projection,...
    'coast',output.coast,'color',output.border_color,'states',output.state,'behrdir',output.behr_dir,...
    'fileprefix',output.behr_prefix,'flags',output.flags,'clouds',output.cloud_type,'cloudfraccrit',output.cloud_max,...
    'rowanomaly',output.rowanomaly);