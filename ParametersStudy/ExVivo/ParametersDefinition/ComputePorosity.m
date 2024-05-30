function [porosity, imageDisplay] = ComputePorosity(binaryImage, boundaryEndost, width, plotBoundary)
% This function allows to compute a percentage of porosity which is
% represented by the volume of bone over the total volume in a specified
% width around the boundary.
% 
% Input :   binaryImage - Binarized image of an ex-vivo X-Rayed image.
%           boundaryEndost - Boundary of the endost, usually extracted with
%           the function ExtractBoundary.
%           nbWavelength - Width of the boundary in which the porosity is
%           computed, in number of wavelength inside the bone.
%           plotBoundary - boolean that indicate weither to plot the
%           extended boundary.
% 
% Output :  porosity - Percentage of bone in the total volume = BV/TV *100
%
% See also : ExpandParabola
    
    % Define the width corresponding to a number of wavelength
    pixelSize = 9e-3; 
    boundaryWidth = width/pixelSize;
    
    [surface, surfaceExtended] = ExpandParabola(boundaryEndost, boundaryWidth);
    
    boneVolume = 0;
    totalVolume = 0;
    imageDisplay = double(binaryImage);

    for i = surfaceExtended(1, surfaceExtended(2,:) < boundaryEndost(2,1))
        if i < boundaryEndost(1,1) || i > boundaryEndost(1,end)
            limSup = max(boundaryEndost(2,:));
        else 
            limSup = max(boundaryEndost(2,boundaryEndost(1,:) == i));
        end
        
        boneVolume = boneVolume + sum(binaryImage(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i));
        totalVolume = totalVolume + numel(binaryImage(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i));

        % Change the value of the pixel in the boundary 
        imageDisplay(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i) =...
            imageDisplay(ceil(surfaceExtended(2,surfaceExtended(1,:) == i)) : limSup, i).*0.5;
    end

    porosity = 100 * (1 - boneVolume/totalVolume);

    if plotBoundary
            figure
            imshow(binaryImage);
            hold on 
            plot(boundaryEndost(1,:), boundaryEndost(2,:), 'r', 'LineWidth', 2)
            hold on         
            plot(surfaceExtended(1,:), surfaceExtended(2,:), 'b', 'LineWidth', 2)
            hold on         
            plot(surface(:,1), surface(:,2), 'g', 'LineWidth', 2)
            xlabel('Width');
            ylabel('Depth');
            grid on
            title('Boundary of the endostium'); 
            legend('Boundary of the endost', 'Limite of the boundary'); hold off;

        fprintf('The bone as a porosity of %.1f %% in a width of %.1f wavelength before the endost thus %.2f mm', ...
            porosity, nbWavelength, nbWavelength*lambda);
    end
end

function [surface, surfaceExtended] = ExpandParabola(boundaryEndost, boundaryWidth)
% This function allows to expand the parabola that fit to the boundary
% endost. 
    surface = FitParabola(boundaryEndost);
    
    % Get coefficient of the parabola
    % In order to compute he model, we reshape it. 
    y = max(surface(:,2)) - surface(:,2);
    x = surface(:,1) - unique(surface(y == max(y), 1));
    coeffParabola = polyfit(x, y, 2); 
    
    % Get maximum and root of the parabola
    r = max(roots(coeffParabola));
    
    % Multiply the coefficients to obtain the extanded parabola 
    A = (coeffParabola(1)*r^2 + boundaryWidth)/(r - boundaryWidth)^2;
    B = coeffParabola(2)*r/(r - boundaryWidth)^2;
    C = coeffParabola(3) - boundaryWidth;
    
    xExtended = (-length(x): 1: length(x))';
    yExtended = A.*xExtended.^2 + B.*xExtended + C;

    % Reshape as previously 
    surfaceExtended(:,1) = xExtended + unique(surface(y == max(y), 1));
    surfaceExtended(:,2) = max(surface(:,2)) - yExtended - 2*boundaryWidth;

    surfaceExtended = surfaceExtended';
end
%% Script to plot the expanded boundary. 
% y = max(surface(:,2)) - surface(:,2);
% x = surface(:,1) - unique(surface(y == max(y), 1));
% angle = atan2(surface(end,2) - surface(1,2), surface(end,1) - surface(1,1));
% 
% % Matrice de rotation
% R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
% 
% % Appliquer la rotation à tous les points de la parabole
% rotated_points = R * [x'; y'];
% 
% % Séparer les points x et y
% x_rotated = rotated_points(1, :);
% y_rotated = rotated_points(2, :);
% 
% x_rotated = x_rotated + unique(surface(y == max(y), 1));
% y_rotated = max(surface(:,2)) - y_rotated;
% 
% x = x + unique(surface(y == max(y), 1));
% y = max(surface(:,2)) - y;
% 
% % Ajuster pour que le début et la fin soient à zéro
% figure;
% plot(x, y, 'b-', 'DisplayName', 'Original Parabola'); % Parabole originale
% hold on;
% plot(x_rotated, y_rotated, 'r-', 'DisplayName', 'Rotated Parabola'); % Parabole pivotée
% legend show;
% xlabel('x');
% ylabel('y');
% title('Rotation de la Parabole');
% grid on;
% 
% figure
% imshow(binaryImage);
% hold on 
% plot(boundaryEndost(1,:), boundaryEndost(2,:), 'r', 'LineWidth', 2)
% hold on         
% plot(surface(:, 1), surface(:,2), 'b', 'LineWidth', 2)
% hold on         
% plot(surface(:, 1), y_rotated)
