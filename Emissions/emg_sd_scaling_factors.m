function [ sd_scales ] = emg_sd_scaling_factors( )
%EMG_SD_SCALING_FACTORS Return the scaling factors for EMG SDs
%   I've gone through and looked at the change in fits when varying the
%   uncertainty for each parameter individually (see notes from June 13,
%   2016) and these are the factors to reduce the uncertainty from
%   minimizing the unexplained variance, rather than sum of squared
%   residuals, that I came up with.  They are returned as the 2x5 matrix:
%
%   [a_lower, x0_lower, mux_lower, sx_lower, B_lower;
%    a_upper, x0_upper, mux_upper, sx_upper, B_upper]
%
%   Josh Laughner <joshlaugh5@gmail.com> 13 June 2016

sd_scales = [0.2, 0.2, 0.4, 0.2, 0.6;
             0.9, 0.6, 0.4, 0.6, 0.2];

end

