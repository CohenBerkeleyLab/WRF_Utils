function [ Var ] = unstagger( Var, dim  )
% U = UNSTAGGER( U, dim ) Converts staggered WRF variables to unstaggered.
%   WRF outputs some variables, notably winds, in staggered coordinates.
%   This is for variables that are calculated on grid edges instead of grid
%   centers. Unstaggering isn't hard, it's just requires a fair bit of
%   typing. This function will unstagger Var along dim by averaging
%   adjacent values along that dimension.
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Mar 2016

E = JLLErrors;
if ~isscalar(dim) || ~isnumeric(dim) || dim > ndims(Var) || dim < 1
    E.badinput('dim must refer to a valid dimension in Var')
end

% We'll put the dimension to unstagger along the first dimension and
% reshape to a 2D matrix to do the unstaggering. Then we'll undo all that
% at the end.
perm_vec = 1:ndims(Var);
perm_vec(dim) = [];
perm_vec = [dim, perm_vec];

Var = permute(Var, perm_vec);
sz = size(Var);
Var = reshape(Var,sz(1),[]);
Var = (Var(1:end-1,:) + Var(2:end,:)) ./ 2;

new_sz = sz;
new_sz(1) = new_sz(1) - 1;

Var = reshape(Var, new_sz);
Var = ipermute(Var, perm_vec);


end

