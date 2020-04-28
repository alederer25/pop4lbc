function [grid, max_dist] = ndgridj(grid_min, grid_max,ns)
% generates a grid and returns all combnations of that grid in a list
% In:
%   grid_min   1  x D  lower bounds of grid for each dimension separately
%   grid_max   1  x D  upper bounds of grid for each dimension separately
%   ns         1  x D  number of points for each dimension separately
% Out:
%   grid       Prod(ns) x D
%   max_dist   1  x 1   maximum distance a point can have in the grid
% Last edited: Armin Lederer, 04/2020

D= numel(ns);

if numel(grid_max) ~= D ||  numel(grid_min) ~= D
    error('grid_max, grid_min and ns must have same dimensions');
end

if ~isvector(grid_max)  ||  ~isvector(grid_min) || ~isvector(ns)
    error('grid_max, grid_min and ns must be vectors');
end

if any(grid_max-grid_min<0)
    error('grid_max ist not always larger than grid_min');
end

gg = cell(1,D);

for i=1:D
    gg{i} = linspace(grid_min(i),grid_max(i),ns(i));
end

if coder.target('MATLAB')
    [gg{1:D}] = ndgrid(gg{1:D});
    
    grid = reshape(permute(cat(D+1, gg{:}),[D:-1:1 D+1]),[],D)';
else
    coder.varsize('grid');
    switch D
        case 1
            ggg1=ndgrid(gg{1});
            grid=ggg1';
        case 2
            [ggg1,ggg2]=ndgrid(gg{1:2});
            grid=reshape(permute(cat(3,ggg1,ggg2),[2:-1:1 3]),[],2)';
        case 3
            [ggg1,ggg2,ggg3]=ndgrid(gg{1:3});
            grid=reshape(permute(cat(4,ggg1,ggg2,ggg3),[3:-1:1 4]),[],3)';
        case 4
            [ggg1,ggg2,ggg3,ggg4]=ndgrid(gg{1:4});
            grid=reshape(permute(cat(5,ggg1,ggg2,ggg3,ggg4),[4:-1:1 5]),[],4)';
        otherwise
            grid=[];
            error('dimension not supported');
    end
end


if nargout >1
    dist = (grid_max-grid_min)./ns;
    max_dist = sqrt(sum(dist.^2))/2;
end


