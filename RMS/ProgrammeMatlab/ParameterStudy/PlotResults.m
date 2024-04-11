clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% PROCESSING OF THE RF DATA 
% Generate a table to stock the geometry and the specularity Map
pathSimuAll = '/calculSSD/salome/Simulation-04avr';
rmsAll = 0.03 + (0:9) * 0.05;      % List of rms values (mm)
corrAll = [0.5 1 2 4];             % List of correlation length (mm)
MapAll = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

%% PLOT ALL GEOMETRY MAP
for i = 1:numel(rmsAll)
    rms = str2double(MapAll.Properties.RowNames{i});
    for j = 1:numel(corrAll)
        % Generate path to simulation directory
        corr = str2double(MapAll.Properties.VariableNames{j});
        format = 'simulation_rms_%.2f_cl_%.1f/';
        simu_dir = [pathSimuAll, sprintf(format, rms, corr)];

        % Recover simulation Map
        [SpecularityMap,~,~] = SimSonic2DReadMap2D([simu_dir, 'Geometry.map2D']);
        MapAll{i, j} = {SpecularityMap};
    end
end

%% PLOT ALL MAPS
figure
t = tiledlayout(numel(rmsAll),numel(corrAll),'TileSpacing','tight');
% Plot Map
for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        % Generate path to simulation directory
        nexttile(t)
        img = MapAll{i,j}{1}(600:1200, :);
        imagesc(img);
%         axis equal
        % Remove x and y axes
        set(gca, 'xtick', [], 'ytick', []);
        
        if i == 1
            title(MapAll.Properties.VariableNames{j}, 'FontSize', 8, 'FontWeight', 'normal');
        end
        if j == 1
            ylabel(MapAll.Properties.RowNames{i}, 'FontSize', 8);
        end
    end
end
title(t, 'Correlation length (mm)', 'FontSize', 10)
ylabel(t, 'Root mean square (mm)', 'FontSize', 10)
xlabel(t, 'Interface profile for various RMS and correlation length', 'FontSize', 15, 'FontWeight', 'bold')

%% PLOT ALL SPECULARITY MAPS 
% Collect all simulations directory for the type of test.
simuDirAll = '/calculSSD/salome/Simulation-04avr'; 

% Create a table to stock the specularity maps
SpecuMapAll = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));
%%
for i = 1:numel(rmsAll)
    rms = str2double(SpecuMapAll.Properties.RowNames{i});
    for j = 1:numel(corrAll)
        % Generate path to simulation directory
        corr = str2double(SpecuMapAll.Properties.VariableNames{j});
        format = 'simulation_rms_%.2f_cl_%.1f';
        simu_dir = fullfile(pathSimuAll, sprintf(format, rms, corr));

        % Recover simulation Map
        if exist(fullfile(simu_dir, 'postProcess.mat'))
            postProcess = load(fullfile(simu_dir, 'postProcess.mat'));
            SpecuMapAll{i, j} = {postProcess.SpecularProbaMap};
        end
    end
end

%% PLOT ALL MAPS
figure
t = tiledlayout(numel(rmsAll),numel(corrAll),'TileSpacing','tight');
% Plot Map
for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        % Generate path to simulation directory
        nexttile(t)
        img = SpecuMapAll{i,j}{1};
        imagesc(img);
%         axis equal
        % Remove x and y axes
        set(gca, 'xtick', [], 'ytick', []);
        
        if i == 1
            title(SpecuMapAll.Properties.VariableNames{j}, 'FontSize', 8, 'FontWeight', 'normal');
        end
        if j == 1
            ylabel(SpecuMapAll.Properties.RowNames{i}, 'FontSize', 8);
        end
    end
end
title(t, 'Correlation length (mm)', 'FontSize', 10)
ylabel(t, 'Root mean square (mm)', 'FontSize', 10)
xlabel(t, 'Specular probability for various RMS and correlation length', 'FontSize', 15, 'FontWeight', 'bold')

%% PLOT ALL THE CORRESPONDING PROBABILITIES 
SpecuProbaAll = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'cell'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

%% COMMON PARAMETERS TO ALL SIMULATION
simuDir = fullfile(pathSimuAll, 'simulation_rms_0.03_cl_0.5');
parameters = load(fullfile(simuDir, 'parameters.mat'));
recorded = LoadRfData(parameters.probe, simuDir);
[acquisition, reconstruction] = GenerateParamRecon(recorded);

%%
for i = 1:numel(rmsAll)
    rms = str2double(SpecuProbaAll.Properties.RowNames{i});
    for j = 1:numel(corrAll)
        % Generate path to simulation directory
        corr = str2double(SpecuProbaAll.Properties.VariableNames{j});
        format = 'simulation_rms_%.2f_cl_%.1f';
        simu_dir = fullfile(pathSimuAll, sprintf(format, rms, corr));

        % Recover simulation Map
        if ~isempty(SpecuMapAll{i, j}{1})
            probability = ProbaROI(SpecuMapAll{i,j}{1}, reconstruction, parameters, false);
            SpecuProbaAll{i,j} = {probability};
        end
    end
end

%% PLOT ALL MAPS
meanROI = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'double'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

stdROI = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'double'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

corrROI = table('Size', [numel(rmsAll), numel(corrAll)], ...
                'VariableType', repmat({'double'}, 1, numel(corrAll)), ...
                'VariableNames', cellstr(string(corrAll)), ...
                'RowNames', cellstr(string(rmsAll)));

figure
t = tiledlayout(numel(rmsAll),numel(corrAll),'TileSpacing','tight');
for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        if ~isempty(SpecuProbaAll{i,j}{1})
            meanROI{i,j} = SpecuProbaAll{i,j}{1}.meanROI;
            stdROI{i,j} = SpecuProbaAll{i,j}{1}.stdROI;
            % corrROI{i,j} = SpecuProbaAll{i,j}{1}.corrROI;

            nexttile(t)
            plot(reconstruction.Xmm, SpecuProbaAll{i,j}{1}.linearROI);
            ylim([0, 1]);
            if i == 1
                title(SpecuMapAll.Properties.VariableNames{j}, 'FontSize', 8, 'FontWeight', 'normal');
            end
            if j == 1
                ylabel(SpecuMapAll.Properties.RowNames{i}, 'FontSize', 8);
            end
        else
            nexttile(t)
        end

    end
end
title(t, 'Correlation length (mm)', 'FontSize', 10)
ylabel(t, 'Root mean square (mm)', 'FontSize', 10)
xlabel(t, 'Specular probability along the lateral position', 'FontSize', 15, 'FontWeight', 'bold')

%%
legendPlot = {};
figure
for i = 1:numel(rmsAll)
    for j = 1:numel(corrAll)
        if ~isempty(SpecuProbaAll{i,j}{1})
            plot(reconstruction.Xmm, SpecuProbaAll{i,j}{1}.linearROI);
            ylim([0, 1]);
            legendPlot{end+1} = ['RMS = ', SpecuProbaAll.Properties.RowNames{i}, ' CORR = ', SpecuProbaAll.Properties.VariableNames{j}];
            
            hold on
        end
    end
end
xlabel('Lateral position (mm)', Interpreter='latex')
ylabel('Specular probability', Interpreter='latex')
title('Specular probability along the lateral position');
legend(legendPlot)
ylim([0, 1]);

 
