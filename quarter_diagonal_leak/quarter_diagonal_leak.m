function quarter_diagonal_leak(h, ht, hl, le, ux, uy, wi, x, y, z)
mrst_startup_path = '/Users/dmar/Github/py-micp/MRST/startup.m';
if exist(mrst_startup_path, 'file') ~= 2
    error('MRST startup file not found: %s', mrst_startup_path);
end
run(mrst_startup_path)
%  Write a grid structure with a diagonal leakage path in the file
%  'GRID.INC' in the 'decks' folder for simulation in
%  OPM Flow
%
% SYNOPSIS:
%       quarter_diagonal_leak(h, ht, hl, le, ux, uy, wi, x, y, z)
%
% PARAMETERS:
%   h    - Scalar cell data (height of the domain).
%
%   ht   - Scalar cell data (height of the top aquifer).
%
%   hl   - Scalar cell data (height of the lower aquifer).
%
%   le    - Scalar cell data (length of the domain).
%
%   ux   - Scalar cell data (x-gap between leak path and left side).
%
%   uy   - Scalar cell data (y-gap between leak path and lower side).
%
%   wi   - Scalar cell data (width of the domain).
%
%   x    - Array cell data (discretization of the x-direction).
%
%   y    - Array cell data (discretization of the y-direction).
%
%   z    - Array cell data (discretization of the z-direction).
%
% NOTES:
%   Function `quarter_diagonal_leak` only creates the diagonal leakage path
%   consisting of one cell aperture.
%
% EXAMPLE:
%   le = 100; h = 30; wi = 10;
%   x = [0 : 21 le * exp(-1.5 : 0.05: 0)];
%   y = [0 : 21 wi * exp(-1.5 : 0.05: 0)];
%   z = [0 : 5 5.5 : 0.5 : h - 5 h - 4 : h];
%   quarter_diagonal_leak(h, 5, 5, le, 18, 11, wi, x, y, z)

    cell_count_x = numel(x) - 1;
    cell_count_y = numel(y) - 1;
    cell_count_z = numel(z) - 1;
    leak_cell_x = sum(x < ux);          % Locate the number of x-cells from the leak
    leak_cell_y = sum(y < uy);          % Locate the number of y-cells from the leak
    top_cell_count = sum(z < ht);       % Locate the number of top cells
    lower_cell_count = sum(z > h - hl); % Locate the number of lower cells
    grdecl = simpleGrdecl([cell_count_x, cell_count_y, cell_count_z], 0, 'flat', true(), 'undisturbed', true(), 'physDims', [2 * le, 2 * wi, h]);
    % Set in ACTNUM as 0 the removed cells
    middle_start = cell_count_x * cell_count_y * top_cell_count + 1;
    middle_stop = cell_count_x * cell_count_y * cell_count_z - cell_count_x * cell_count_y * lower_cell_count;
    if middle_start <= middle_stop
        grdecl.ACTNUM(middle_start : middle_stop, 1) = 0;
    end
    leak_layers = top_cell_count : cell_count_z - lower_cell_count;
    diagonal_cell_count = leak_cell_x - leak_cell_y;
    if diagonal_cell_count > 0
        diagonal_steps = 0 : diagonal_cell_count - 1;
        diagonal_row_offsets = (diagonal_steps + leak_cell_y) * cell_count_x;
        diagonal_col_offsets = leak_cell_x - diagonal_steps;
        first_diagonal_indices = diagonal_row_offsets + diagonal_col_offsets;
        second_diagonal_indices = first_diagonal_indices + 1;
        base_diagonal_indices = unique([first_diagonal_indices, second_diagonal_indices]);
        layer_offsets = cell_count_x * cell_count_y * leak_layers;
        leak_indices = bsxfun(@plus, layer_offsets(:), base_diagonal_indices(:).');
        grdecl.ACTNUM(unique(leak_indices(:)), 1) = 1;
    end
    % Create the GRID.INC file
    grdecl2 = grdecl;
    node_count_x = cell_count_x + 1;
    node_count_y = cell_count_y + 1;
    pillar_count = node_count_x * node_count_y;
    x_coords = repmat(x(:), node_count_y, 1);
    y_coords = kron(y(:), ones(node_count_x, 1));
    grdecl2.COORD = zeros(6 * pillar_count, 1);
    grdecl2.COORD(1 : 6 : end) = x_coords;
    grdecl2.COORD(2 : 6 : end) = y_coords;
    grdecl2.COORD(3 : 6 : end) = z(1);
    grdecl2.COORD(4 : 6 : end) = x_coords;
    grdecl2.COORD(5 : 6 : end) = y_coords;
    grdecl2.COORD(6 : 6 : end) = z(end);
    zcorn_face_size = cell_count_x * cell_count_y * 4;
    zcorn_layer_size = cell_count_x * cell_count_y * 8;
    grdecl2.ZCORN = zeros(size(grdecl2.ZCORN));
    grdecl2.ZCORN(1 : zcorn_face_size) = z(1);
    if cell_count_z > 1
        middle_zcorn = reshape(repmat(z(2 : cell_count_z), zcorn_layer_size, 1), [], 1);
        grdecl2.ZCORN(zcorn_face_size + 1 : end - zcorn_face_size) = middle_zcorn;
    end
    grdecl2.ZCORN(end - zcorn_face_size + 1 : end) = z(end);
    writeGRDECL(grdecl2, 'decks/GRID.INC');
end

%{
Copyright 2021-2026, NORCE Research AS, Computational
Geosciences and Modeling.
This file is part of the py-micp module.
py-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
py-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
