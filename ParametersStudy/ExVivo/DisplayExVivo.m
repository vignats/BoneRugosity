% COMPUTE PARAMETERS OF EX-VIVO BONES
clear; 
close all;
clc;
addpath(genpath('~/Documents')); 
addpath(genpath('/calculSSD/salome')); 

%% 
bones = {'245D', '227G', '267G'};
slices = {{1134, 3195, 3852, 5511}, {2002, 3595, 5614, 6721}, {1700, 3410, 5478, 6716}}; 
% slicesbis = {{1134, 3195, 3852, 5511}, {1590, 3595, 5614, 6721}, {1700, 3410, 5478, 6716}}; 
kcAll = [0.041 0.046 0.057 0.064; 0.051 0.046 0.033 0.037; 0.033 0.027 0.030 0.027];
figure
t = tiledlayout(numel(bones), numel(slices{1}), 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:numel(bones)  
    bone.id = bones{i};         % Bone from ex-vivo files
    for j = 1:numel(slices{i})
        bone.image = slices{i}{j};        % Slice selected
        % % imagebis = slicesbis{i}{j}; 
        % 
        % Load data 
        simulation_name = ['Bone', bone.id, '-Image', sprintf('%04d', bone.image), '/'];
        load(fullfile('/calculSSD/salome/Simulation-10mai/', simulation_name, 'parameters.mat'));

        % Binarize bone image
        dirname = '/calculSSD/salome/BoneImage/';
        % dirname = ['/calculSSD/Dossier partagé image os exvivo/', bone.id, '/ZONE_US_0', int2str(j)];
        file = ['SAMPLE_', bone.id, '_SLICE_', sprintf('%04d', bone.image), '.bmp']; 
        % file = ['SAMPLE_', bone.id, '_SLICE_', sprintf('%04d', imagebis), '.bmp']; 
        filename = fullfile(dirname, file); 

        bone_bmp = imread(filename); 
        threshold = graythresh(bone_bmp); % Find an automatic threshold
        binaryImage = imbinarize(bone_bmp, threshold);

        if i == 1 && j ~=1
            if j == 2 || j == 4
                angle = -10;
            elseif j == 3
                angle = -15;
            end
            binaryImage = imrotate(binaryImage, angle,'bilinear','loose');
            binaryImage = binaryImage(1:find(...
            sum(binaryImage(: , 1: size(binaryImage, 2)/2),2) > 10, 1, 'last'), :);
        end 
        % 
        % % Extract endost
        lambda = (medium.cp(2) / (probe.fc))*1e3; 
        % kc = 1/(lambda*50); % Need to be change to correspond to a wavelength in the bone with the appropriate velocity. 
        [profile, roughness, ~, kc] = FilterProfil(binaryImage, kcAll(i,j));
        % fprintf(['\nThe cut off frequency for the bone ' bone.id ' and the slices %04d is %.03f mm -1'], bone.image,kc)

        % Compute the correlation length and rms of the roughness profile 
        Rms = rms(roughness);
        corr = ComputeCorr(roughness, grid.step);

        % Compute the porosity in a surface of one wavelength around the boundary
        nbWavelength = 1/2;
        width = lambda * nbWavelength;
        [boundaryEndost, ~] = ExtractBoundary(binaryImage);
        [porosity, imageDisplay] = ComputePorosity(binaryImage, boundaryEndost, width, false);

        % Plot results
        nexttile(t)
        stepXRay = 0.009; % Step of the X-Ray image (mm)
        X = 0:stepXRay:size(imageDisplay, 2)*stepXRay; %X = X - mean(X);
        Z = flipud(0:stepXRay:size(imageDisplay, 1)*stepXRay);
        imagesc(X, Z, imageDisplay)
        hold on
        plot(boundaryEndost(1,:)*stepXRay, boundaryEndost(2,:)*stepXRay, 'r', 'LineWidth', 2)
        xlabel('Width (mm)');
        ylabel('Depth (mm)');
        hold off
        axis image

        % Calculate dimensions of current subplot
        ax = gca;
        set(ax, 'Units', 'normalized'); 
        pos = ax.Position;

        % Adjust position of annotation box
        dim = [pos(1), pos(2) + pos(4) - 0.08, 0.3, 0.1];
        str = {[texlabel('sigma') sprintf(' = %.2fmm', Rms)]...
            [texlabel('rho') sprintf(' = %.3fmm', corr)] sprintf('E.Por = %.1f%%', porosity)};
        annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor','w',...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

    end
end
title(t, 'Parameters of Bone Health on Ex-Vivo X-Ray Images', 'FontSize', 10)
