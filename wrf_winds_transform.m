function [ Ue, Ve ] = wrf_winds_transform( U, V, COSALPHA, SINALPHA )
% WRF_WINDS_TRANSFORM Transform winds to earth-relative orientation
%   Function that will correctly calculate U and V wind components of
%   winds given the COSALPHA and SINALPHA variables from WRF output,
%   handling the conversion from grid-relative winds to
%   earth-relative winds.  Normally there is very little difference
%   between the two.
%
%   Note that this will assume that COSALPHA and SINALPHA will be the
%   same for each 2D slice of U and V and so will only use the first
%   2D slice of COSALPHA and SINALPHA. It will also unstagger U and
%   V.

sz_U = size(U);
sz_V = size(V);
stag_U = zeros(size(sz_U));
stag_U(1) = 1;
stag_V = zeros(size(sz_V));
stag_V(2) = 1;

if ndims(U) ~= ndims(V) || (~all(size(U) == size(V)) && ~all(size(U)-stag_U == size(V) - stag_V))
    E.badinput('U and V should be the same size (they may be left staggered).')
elseif ndims(COSALPHA) ~= ndims(SINALPHA) || ~all(size(COSALPHA) == size(SINALPHA))
    E.badinput('COSALPHA and SINALPHA should be the same size')
end

% Unstagger if needed
if size(U,1)-1 == size(V,1)
    U = unstagger(U,1);
end
if size(V,2)-1 == size(U,2)
    V = unstagger(V,2);
end

% Warn if the alphas have dimensions that won't be used
if ndims(COSALPHA) > 2 || ndims(SINALPHA) > 2 %#ok<*ISMAT>
    warning('Only the first 2D slice of COSALPHA and SINALPHA will be used - this is assumed to be unchanged for all the WRF output')
    COSALPHA = COSALPHA(:,:,1);
    SINALPHA = SINALPHA(:,:,1);
end

% Loop through and apply the transformation to each 2D slice of the
% wind
Ue = nan(size(U));
Ve = nan(size(V));
for a=1:prod(sz_U(3:end))
    Ue(:,:,a) = U(:,:,a) .* COSALPHA - V(:,:,a) .* SINALPHA;
    Ve(:,:,a) = V(:,:,a) .* COSALPHA + U(:,:,a) .* SINALPHA;
end

end

function M = unstagger(M, dim)
permvec = 1:ndims(M);
permvec(1) = dim;
permvec(dim) = 1;
M = permute(M,permvec);
M = (M(1:end-1,:,:,:,:) + M(2:end,:,:,:,:))/2;
M = permute(M, permvec);
end