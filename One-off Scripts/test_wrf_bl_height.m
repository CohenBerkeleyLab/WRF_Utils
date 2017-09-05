function [ fraction_good, fraction_ok, fraction_undefined, indices_good, indices_ok, indices_bad, indices_undefined  ] = test_wrf_bl_height( z, no2, is_z_pres, pblh )
%TEST_WRF_BL_HEIGHT Plots tests of ways to find the chemical BL height from WRF profiles
%   Detailed explanation goes here


E = JLLErrors;

if ndims(z) ~= 3 || ndims(no2) ~= 3
    E.badinput('Both input matrices must be 3D')
elseif any(size(z) ~= size(no2)) && size(z,3) ~= size(no2,3)+1
    E.badinput('Both input matrices must have the same size, or z must have at most one more element in the third dimension')
end

if ~exist('is_z_pres','var')
    is_z_pres = false;
end
if ~exist('pblh','var')
    use_wrf_blh = false;
else
    use_wrf_blh = true;
end

if size(z,3) == size(no2,3)+1
    z = z(:,:,1:end-1);
end

% This function will randomly pick profiles from the matrix provided,
% calculate the boundary layer height, and plot both the profile and the
% height. It will then ask if the height is accurate and keep track of what
% percent (and which ones) are good or bad. This will continue until the
% user says to quit.

n = 0;
n_good = 0;
n_ok = 0;
n_undefined = 0;
indices_good = {};
indices_ok = {};
indices_bad = {};
indices_undefined = {};


while true
    n = n+1;
    i = randi(size(z,1));
    j = randi(size(z,2));
    
    z_slice = squeeze(z(i,j,:));
    no2_slice = squeeze(no2(i,j,:));
    
    f=figure;
    line(no2_slice, z_slice, 'linewidth', 2, 'color', 'k', 'linestyle', '-');
    if is_z_pres
        set(gca,'ydir','reverse');
    end
    
    method = 'exp2';
    
    while true
        if ~use_wrf_blh
            blh = find_bdy_layer_height(no2_slice, z_slice, method, 'altispres', is_z_pres);
            l=line([min(no2_slice(:)), max(no2_slice(:))], [blh blh], 'linewidth', 2, 'color', 'b', 'linestyle', '--');
            d = (no2_slice(1) - median(no2_slice))/no2_slice(1) * 100;
            title(sprintf('%.2f',d));
        else
            l = line([min(no2_slice(:)), max(no2_slice(:))], [pblh(i,j)+z_slice(1), pblh(i,j)+z_slice(1)]);
        end
        
        s = ask_question;
        
        switch s
            case 'y'
                n_good = n_good + 1;
                indices_good = cat(1,indices_good, {[i,j],method});
                break
            case 'n'
                indices_bad = cat(1, indices_bad, {[i,j],method});
                break
            case 'o'
                n_ok = n_ok + 1;
                indices_ok = cat(1,indices_ok, {[i,j],method});
                break
            case 'e'
                if ~use_wrf_blh && ~strcmp(method,'exp')
                    delete(l);
                    if strcmp(method, 'exp2')
                        method = 'exp1.5';
                    elseif strcmp(method, 'exp1.5')
                        method = 'exp';
                    end
                else
                    n_undefined = n_undefined + 1;
                    indices_undefined = cat(1, indices_undefined, {[i,j],method});
                    break
                end
            case 'q'
                n = n-1; % we never categorized this last profile.
                break
        end
    end
    
    close(f);
    
    if strcmp(s,'q')
        break
    end
end

fraction_good = n_good / n;
fraction_ok = n_ok / n;
fraction_undefined = n_undefined / n;

end

function s = ask_question
fprintf('Is the BL height for this profile accurate?\n')
quest = '   Enter ''y'' for yes, ''n'' for no, ''o'' for okay, ''e'' if no BL height shown, or ''q'', to quit: ';
while true
    s = lower(input(quest, 's'));
    s = s(1);
    if ~ismember(s,'ynoeq')
        fprintf(' * Please enter y, n, o, e, or q! *\n')
    else
        break
    end
end
end

