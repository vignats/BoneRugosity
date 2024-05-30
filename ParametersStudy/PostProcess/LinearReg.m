clear; 
close all;
clc;
addpath(genpath('~/Documents'));
addpath(genpath('/calculSSD/salome'));

%% Formated data in vector
simulationDate = {'19avr', '24avr', '29avr'};

% Creation of input matrice, the column corresponds to 
% Correlation length    Root Mean Square    Porosity    Diameter of pores (mm)
inputAll = zeros([length(simulationDate)*40, 2]);
outputAll = zeros([length(simulationDate)*40, 1]);

rmsAll = 0.03 + (0:9) * 0.05;      % List of rms values (mm)
corrAll = [0.5 1 2 4];             % List of correlation length (mm)

inputCell = cell([length(simulationDate), 1]);
outputCell = cell([length(simulationDate), 1]);
ResultsAll = cell([length(simulationDate), 1]);
for simu = 1:length(simulationDate)
    load(fullfile('/calculSSD/salome', ['Simulation-' simulationDate{simu}], 'Results.mat'));
    inputSpecular = zeros([40, 2]);     % 1st column is correlation length, 2nd is rms
    outputSpecular = zeros([40, 1]);
    
    ResultsAll{simu} = Results;

    output = table2array(Results.Stat.meanROI);
    
    for i = 1:length(corrAll)
        idx = (i-1)*10;
        for j = 1:length(rmsAll)
            InitParam = Results.InitParam;
            inputSpecular(idx + j, 1) = InitParam.Corr{j,i};
            inputSpecular(idx + j, 2) = InitParam.Rms{j,i};
            
            outputSpecular(idx + j, 1) = output(j,i);
        end
    end
    inputAll(40*(simu-1) + 1 : 40*simu, 1:2) = inputSpecular;
    outputAll(40*(simu-1) + 1 : 40*simu) = outputSpecular;
    
    % Save data to plot surfaces in function of real input value
    inputCell{simu} = inputSpecular;
    outputCell{simu} = outputSpecular;
end
%% Linear Regression with both parameters
mdl = fitlm(inputAll,outputAll);

%% Correlation coefficient
coeffAll = corrcoef([outputAll inputAll]);

% With the ratio of correlation over rms
ratioAll = inputAll(:,2)./ inputAll(:,1);
coeffRatioAll = corrcoef([outputAll ratioAll]);

%% Correaltion coefficient for each data 
coeffOneLayer = corrcoef([outputAll(1:40) inputAll(1:40, :)]);
coeffTwoLayer = corrcoef([outputAll(41:80) inputAll(41:80, :)]);
coeffTwoLayer10pores = corrcoef([outputAll(81:120) inputAll(81:120, :)]);
%% Plot mean ROI
figure;
% surf(corrAll, rmsAll, table2array(ResultsAll{1}.Stat.meanROI), 'FaceColor', 'red', 'FaceAlpha', 0.5);
% hold on
% surf(corrAll, rmsAll, table2array(ResultsAll{2}.Stat.meanROI), 'FaceColor', 'blue', 'FaceAlpha', 0.5);
% hold on
% surf(corrAll, rmsAll, table2array(ResultsAll{3}.Stat.meanROI), 'FaceColor', 'green', 'FaceAlpha', 0.5);
% hold on 
scatter3(inputCell{1}(:,1), inputCell{1}(:,2),outputCell{1}, 'r*');
hold on;
scatter3(inputCell{2}(:,1), inputCell{2}(:,2), outputCell{2}, 'b*');
hold on;
scatter3(inputCell{3}(:,1), inputCell{3}(:,2), outputCell{3}, 'g*');

xlabel('Correlation length (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability');
title('Mean Specular Probability in the Region of Interest');
legend('No pores', '10um pores', '30um pores');
%% Plot with real inputs
% Augmenter la résolution de la grille
numPoints = 40; % Par exemple, 50 points dans chaque direction

[X1, Y1] = meshgrid(linspace(min(inputCell{1}(:,1)), max(inputCell{1}(:,1)), numPoints),...
    linspace(min(inputCell{1}(:,2)), max(inputCell{1}(:,2)), numPoints));
Z1 = griddata(inputCell{1}(:,1), inputCell{1}(:,2), outputCell{1}, X1, Y1, 'natural');

[X2, Y2] = meshgrid(linspace(min(inputCell{2}(:,1)), max(inputCell{2}(:,1)), numPoints),...
    linspace(min(inputCell{2}(:,2)), max(inputCell{2}(:,2)), numPoints));
Z2 = griddata(inputCell{2}(:,1), inputCell{2}(:,2), outputCell{2}, X2, Y2, 'natural');

[X3, Y3] = meshgrid(linspace(min(inputCell{3}(:,1)), max(inputCell{3}(:,1)), numPoints),...
    linspace(min(inputCell{3}(:,2)), max(inputCell{3}(:,2)), numPoints));
Z3 = griddata(inputCell{3}(:,1), inputCell{3}(:,2), outputCell{3}, X3, Y3, 'natural');

figure;
surf(X1, Y1, Z1, 'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
surf(X2, Y2, Z2, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
surf(X3, Y3, Z3, 'FaceColor', 'green', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
scatter3(inputCell{1}(:,1), inputCell{1}(:,2), outputCell{1}, 'r*');
hold on;
scatter3(inputCell{2}(:,1), inputCell{2}(:,2), outputCell{2}, 'b*');
hold on;
scatter3(inputCell{3}(:,1), inputCell{3}(:,2), outputCell{3}, 'g*');

%% Plot the scatter map
% Extract numeric data from cell arrays
x = inputCell{1}(:,1); % First column of inputs
y = inputCell{1}(:,2); % Second column of inputs (if exists)
probaRoi = outputCell{1}; % Output values

% Fit the data using a polynomial surface fit
sf = fit([x, y], probaRoi, "linearinterp",ExtrapolationMethod="linear");

% Create a 3D scatter plot with the fitted surface
figure;
plot(sf, [x, y], probaRoi); % Scatter plot of original data
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Scatter Plot with Fitted Surface');
grid on;
%%
figure;
scatter3(inputCell{1}(:,1), inputCell{1}(:,2),outputCell{1}, 'red');
hold on;
scatter3(inputCell{2}(:,1), inputCell{2}(:,2), outputCell{2}, 'blue');
hold on;
scatter3(inputCell{3}(:,1), inputCell{3}(:,2), outputCell{3}, 'green');
xlabel('Correlation length (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability');
title('Mean Specular Probability in the Region of Interest');
legend('No pores', '10um pores', '30um pores');

%% Plot in 2D
figure;
col = {{'#E57373', '#F44336', '#D32F2F', '#B71C1C'}, {'#7986CB  ', '#3F51B5', '#303F9F', '#1A237E'},...
    {'#AED581', '#8BC34A', '#689F38', '#33691E'}};
obj = {'*', 'o', 'd', 'x'};
simulParam = {'0', '10', '30'};

handlesColor = [];
labelsColor = {};
handlesShape = [];
labelsShape = {};

for i = 1:length(ResultsAll)
    for j = 1:length(obj)
        s = scatter(table2array(ResultsAll{i}.InitParam.Rms(:,j)), table2array(ResultsAll{i}.Stat.meanROI(:,j)),...
            obj{j}, 'MarkerEdgeColor', col{i}{j});

        [curve, ~] = fit(table2array(ResultsAll{i}.InitParam.Rms(:,j)), table2array(ResultsAll{i}.Stat.meanROI(:,j)),'poly1');
        p = plot(curve);
        p.Color = col{i}{j};
        % p.LineWidth = 2;
        
        if i == 1
            handlesShape(end+1) = s;
            labelsShape{end+1} = sprintf('%s = %.1f mm', texlabel('rho'), corrAll(j));
        end
        
        if j ==1
            handlesColor(end+1) = p;
            labelsColor{end+1} = sprintf('Por.Size = %s %sm', simulParam{i}, texlabel('mu'));
        end
        hold on;
    end
end
hold off
handles = [handlesColor, handlesShape];
labels = [labelsColor, labelsShape];

xlabel('Root mean square (mm)');
ylabel('Mean Specular Probability');
title('Mean Specular Probability in the Region of Interest');
% legend(handlesColor, labelsColor, 'Location', 'northeastoutside');
% legend(handlesShape, labelsShape, 'Location', 'northeastoutside');
%% Plot the mean value for rms
figure;
plot(table2array(mean(ResultsAll{1}.InitParam.Corr, 2)), table2array(mean(ResultsAll{1}.Stat.meanROI, 2)), 'r*');
hold on;
plot(table2array(mean(ResultsAll{1}.InitParam.Corr, 2)), table2array(mean(ResultsAll{2}.Stat.meanROI, 2)), 'b*');
hold on;
plot(table2array(mean(ResultsAll{1}.InitParam.Corr, 2)), table2array(mean(ResultsAll{3}.Stat.meanROI, 2)), 'g*');
xlabel('Root mean square (mm)');
ylabel('Mean Specular Probability');
title('Mean Specular Probability in the Region of Interest');
legend('No pores', '10um pores', '30um pores');

%%
figure;
hold on;

% Couleurs pour chaque groupe
colors = {'red', 'blue', 'green'};
labels = {'No pores', '10um pores', '30um pores'};

% Boucle sur chaque groupe de points
for i = 1:length(inputCell)
    corrValues = inputCell{i}(:,1);
    rmsValues = inputCell{i}(:,2);
    probaRoi = outputCell{i};

    % Tracer le nuage de points avec des croix
    scatter3(corrValues, rmsValues, probaRoi, 'x', 'MarkerEdgeColor', colors{i});

    % Ajuster le plan pour le groupe
    X = [ones(size([corrValues, rmsValues], 1), 1), [corrValues, rmsValues]];
    b = X \ probaRoi;

    % Tracer le plan pour le groupe
    [xGrid, yGrid] = meshgrid(linspace(min(corrValues), max(corrValues), 10), linspace(min(rmsValues), max(rmsValues), 10));
    zGrid = b(1) + b(2) * xGrid + b(3) * yGrid;
    surf(xGrid, yGrid, zGrid, 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', colors{i});
end

% Labels and title
xlabel('Correlation length (mm)');
ylabel('Root mean square (mm)');
zlabel('Mean Specular Probability');
title('Mean Specular Probability in the Region of Interest');

% Créer une légende en utilisant les handles des scatter et surf
legend('E.Pore = 0', 'Fitted surface, E.Pore = 0', 'E.Pore = 10mm', 'Fitted surface, E.Pore = 10mm',...
    'E.Pore = 30mm', 'Fitted surface, E.Pore = 30mm');

hold off;


