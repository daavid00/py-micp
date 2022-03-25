function quarter_double_leak(H, ht, hl, L, ux1, ux2, uy1, uy2, Wi, X, Y, Z)
run '~/Github/mrst-2021b/startup.m'
%  Write a grid structure with two straight leakage paths in the file
%  'quarter_double_leak.grdecl' in the 'inputs' folder for simulation in
%  OPM Flow
%
% SYNOPSIS:
%       quarter_double_leak(H, ht, hl, L, ux1, ux2, uy1, uy2, Wi, X, Y, Z)
%
% PARAMETERS:
%   H    - Scalar cell data (height of the domain).
%
%   ht   - Scalar cell data (height of the top aquifer).
%
%   hl   - Scalar cell data (height of the lower aquifer).
%
%   L    - Scalar cell data (length of the domain).
%
%   ux1  - Scalar cell data (X-gap between leak path and left side).
%
%   ux2  - Scalar cell data (X-gap between leak path and left side).
%
%   uy1  - Scalar cell data (Y-gap between leak path and lower side).
%
%   uy2  - Scalar cell data (Y-gap between leak path and lower side).
%
%   Wi   - Scalar cell data (width of the domain).
%
%   X    - Array cell data (discretization of the x-direction).
%
%   Y    - Array cell data (discretization of the y-direction).
%
%   Z    - Array cell data (discretization of the z-direction).
%
% NOTES:
%   Function `quarter_double_leak` only creates two straigth leakage paths
%   consisting of one cell aperture.
%
% EXAMPLE:
%   L = 100; H = 30; Wi = 10;
%   X = [0 : 21 L * exp(-1.5 : 0.05: 0)];
%   Y = [0 : 21 Wi * exp(-1.5 : 0.05: 0)];
%   Z = [0 : 5 5.5 : 0.5 : H - 5 H - 4 : H];
%   quarter_double_leak(H, 5, 5, L, 11, 18, 11, 18, Wi, X, Y, Z)

    xx = size(X, 2) - 1;
    yy = size(Y, 2) - 1;
    zz = size(Z, 2) - 1;
    G = tensorGrid(X, Y, Z);
    GcD1 = G.cartDims(1);
    GcD2 = G.cartDims(2);
    GcD3 = G.cartDims(3);
    Gnc = G.nodes.coords;
    ix1 = sum(X < ux1);        % Locate the number of x-cells from the leak
    ix2 = sum(X < ux2);        % Locate the number of x-cells from the leak
    iy1 = sum(Y < uy1);        % Locate the number of y-cells from the leak
    iy2 = sum(Y < uy2);        % Locate the number of y-cells from the leak
    it = sum(Z < ht);          % Locate the number of top cells
    il = sum(Z > H - hl);      % Locate the number of lower cells
    grdecl = simpleGrdecl([xx, yy, zz], 0, 'flat', true(), ...
                    'undisturbed', true(), 'physDims', [2 * L, 2 * Wi, H]);
    % Set in ACTNUM as 0 the removed cells
    for i = xx * yy * it + 1 : xx * yy * zz - xx *yy * il
        grdecl.ACTNUM(i, 1) = 0;
    end
    for k = it : zz - il
        grdecl.ACTNUM(xx * yy * k + iy1 * xx + ix1, 1) = 1;
        grdecl.ACTNUM(xx * yy * k + iy2 * xx + ix2, 1) = 1;
    end
    % Create the .GRDECL file
    grdecl2 = grdecl;
    j = 1;
    for k = 0 : GcD2
        for i = 1 : GcD1 + 1
            grdecl2.COORD(j) = Gnc(i, 1);
            grdecl2.COORD(j + 1) = Gnc(k * (GcD1 + 1) +1, 2);
            grdecl2.COORD(j + 2) = Gnc(1, 3);
            grdecl2.COORD(j + 3) = Gnc(i, 1);
            grdecl2.COORD(j + 4) = Gnc(k * (GcD1 + 1) + 1, 2);
            grdecl2.COORD(j + 5) = Gnc(end, 3);
            j = j + 6;
        end
    end
    grdecl2.ZCORN(1 : GcD1 * GcD2 * 4) = Gnc(1, 3);
    for i = 1 : GcD3 - 1
        grdecl2.ZCORN(GcD1 * GcD2 * 4 + 1 + GcD1 * GcD2 * 8 * (i - 1) : ...
                     GcD1 * GcD2 * 8 * i + GcD1 * GcD2 * 4) = Gnc((i) * ...
                                     ((GcD1 + 1) * (GcD2 + 1) + 1) + 1, 3);
    end
    grdecl2.ZCORN(end - GcD1 * GcD2 * 4 + 1 : end) = Gnc(end, 3);
    writeGRDECL(grdecl2, 'inputs/quarter_double_leak.grdecl');
end
%{
Copyright 2021-2022, NORCE Norwegian Research Centre AS, Computational
Geosciences and Modeling.
This file is part of the py-micp module.
py-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
