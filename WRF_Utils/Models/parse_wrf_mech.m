function [ J, species, isfixed ] = parse_wrf_mech(  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

species_file = '/Users/Josh/Documents/MATLAB/BEHR/WRF_Utils/Models/R2SMH/r2smh.spc';
eqn_file = '/Users/Josh/Documents/MATLAB/BEHR/WRF_Utils/Models/R2SMH/r2smh.eqn';

% The mechanism will be constructed as a vector of anonymous functions that
% accept a concentration vector and photolysis rate vector.
%
% An additional emissions vector will be added at each timestep as well.

[species, isfixed] = parse_species(species_file);
J = construct_mechanism(eqn_file, species);
end

function [species, isfixed] = parse_species(species_file)
% This function will return the list of species defined for the mechanism
% as well as a logical vector indicating if those species are to be fixed
% by some physical process rather than the chemistry.
%
% First read in the species names from the species file, as well as if it
% is meant to be solved by chemistry or fixed by physical processes.
species = cell(1,1000); % will cut down in the end
isfixed = false(1,1000);
fixed_i = 0;
i=1;
fid = fopen(species_file);
tline = fgetl(fid);
while ischar(tline)
    tline = strtrim(tline);
    if ismember('#DEFVAR',tline) % the following species are variable
        fixed_i = 0;
    elseif ismember('#DEFFIX',tline) % the following species are fixed
        fixed_i = 1;
    elseif ismember('#',tline) % unanticipated definition bloc
        fprintf('Unknown definition block: %s\n',tline);
    else
        % Some lines have curly braces. Since I don't know exactly what
        % they mean, and I'm only interested in what species are needed,
        % they can be removed.
        parsed_line = strrep(tline,'{','');
        parsed_line = strrep(parsed_line,'}','');
        parsed_line = strsplit(parsed_line,'='); % species comes before the = sign
        spec = strtrim(parsed_line{1});
        xx = strcmp(spec, species);
        if sum(xx) > 0 && isfixed(xx) ~= fixed_i
            error('species_definition:species_fixed_and_var','Specie %s is defined as both a fixed and variable specie',spec)
        elseif sum(xx) > 0 
            fprintf('Redundant definition of %s, skipping\n',spec);
        else
            species{i} = spec;
            isfixed(i) = fixed_i;
            i=i+1;
        end
    end
    tline = fgetl(fid);
end
species(i:end) = [];
isfixed(i:end) = [];
end

function J = construct_mechanism(eqn_file, species)
J = cellmat(numel(species),1,1,1);

fid = fopen(eqn_file);
tline = fgetl(fid);
while ischar(tline)
    % skip comment lines or "dummy" reactions 
    if ismember('#',tline) || ~isempty(regexp(tline,'JUNK','once'))
    else
        % Get rid of the equation number which is enclosed in {}
        i=strfind(tline,'}');
        tline=tline(i+1:end);
        
        % Split into equation and rate constant
        tmp = strtrim(strsplit(tline,':'));
        rate_const = tmp{2};
        % then products and reactants. First remove additional curly braces
        % which seem to be used to define fixed species, but we already
        % know that from the species list.
        tmp{1} = regexprep(tmp{1},'[{}]','');
        tmp = strtrim(strsplit(tmp{1},'='));
        products = strtrim(strsplit(tmp{2},'+'));
        reactants = strtrim(strsplit(tmp{1},'+'));
        
        % Handle coefficients in reactants and products, or a product being
        % specified multiple times.
        
        % Parse the rate constant
        k=1; % temp debugging code.
        
        % Construct the derivative function and add it to any existing such
        % functions in the mechanism for this species. Remove 'hv' from
        % the reactants at this point (photons)
        reactants = reactants(~strcmp(reactants,'hv'));
        rr = ismember(species, reactants);
        if sum(rr) ~= numel(reactants)
            nf = ~ismember(reactants, species);
            error('mechanism_parse:unknown_reactant','The reactant %s cannot be identified in the species list',strjoin(reactants(nf),', '));
        end
        
        f = sprintf('%e*prod(c(%s))', k, mat2str(find(rr)));
        
        % Add this for each product
        pp = find(ismember(species, products));
        if numel(pp) ~= numel(products)
            nf = ~ismember(products, species);
            error('mechanism_parse:unknown_reactant','The product %s cannot be identified in the species list',strjoin(products(nf),', '));
        end
        for j=pp
            if isa(J{j},'function_handle')
                J{j} = eval(sprintf('@(c,j) %s + %s',f,func2str(J{j})));
            else
                J{j} = eval(sprintf('@(c,j) %s',f));
            end
        end
    end
    tline = fgetl(fid);
end
end