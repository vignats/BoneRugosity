% Script to visualize the snapshot of a simulation 
% Parameters computation
simulationPath = '/calculSSD/salome/Simulation-pres/simulation_rms_0.00_cl_0.0';
tx = 1; 

to_px = @(mm) round(mm/grid.step); 
[Map,Nx,Nz] = SimSonic2DReadMap2D(fullfile(simulationPath, 'Geometry.map2D'));

X = 0:grid.step:grid.width-grid.step; X = X - mean(X);
Z = flipud(0:grid.step:grid.depth-grid.step);
snapNb = param.length/0.5;      % Number of snap per transmitter, one stap recorded every 0.5 us. . 

% profile = zeros(1, size(Map, 2));
% for i = 1:size(Map, 2)
%     profile(1, i) = find(Map(:, i) == 1, 1, 'last'); % Find last bone element that correspond to the bone ST interface
% end
%% Plot map
figure;
for i = 1:snapNb
    snapName = sprintf('tx_%02d/T11_%03d.snp2D',tx, i);
    snap = SimSonic2DReadSnp2D(fullfile(simulationPath, snapName));
    
    % Plot
    imagesc(X, Z, snap.Data)
    % hold on 
    % plot(X, profile*grid.step, 'w', 'LineWidth', 2)
    hold on
    axis ij image, %shading interp
    pause(0.1) % Pause between the display of each snap
end
