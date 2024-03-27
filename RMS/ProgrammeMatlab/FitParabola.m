function [surface] = FitParabola(boundary_endost)
% FITPARABOLA Fits a parabolic curve to the provided boundary points of the endosteum.
%
%   surface = FitParabola(binary_image, boundary_endost_init) fits a parabolic curve to
%   the boundary points of the endosteum obtained from the binary image. It returns
%   the coordinates of the fitted parabola.
%
%   Inputs:
%   - binary_image: Binary image containing the bone structure.
%   - boundary_endost: Coordinate of the endosteum boundary.
%
%   Output:
%   - surface: Matrix containing the coordinates of the fitted parabola.
%
%   Example:
%       surface = FitParabola(binary_image, boundary_endost)
%
%   See also: fit

    model = fit(boundary_endost(:, 1),  boundary_endost(:, 2), 'poly2');
    y = feval(model, boundary_endost(:, 1));  % Y coordinate of the parabola
    surface = [boundary_endost(:, 1), y];     % Associate the Y coordinate of the parabola with the X coordinate
end