function [max_z_index,max_x_index] = segmentation_Dijkstra_pour_Cecile(im,vertical_freedom_pts)

vertical_freedom_pts = vertical_freedom_pts - 1;

% figure(5577)
% imagesc(im)
% hold on
% xlabel('azimutal distance (pixel)')
% ylabel('range distance/depth (pixel)')
% colorbar
% title('draw polygon for segmentation')
% drawnow
% 
% h = impoly(gca);
% wait(h);
% mask = createMask(h);
% im(~mask) = 0;
% 
% figure(5577)
% imagesc(im)
% hold on
% xlabel('azimutal distance (pixel)')
% ylabel('range distance/depth (pixel)')
% colorbar
    % title('select rectangular region you want to mute')
% drawnow
    
    % MeritFun will contain the final cost for a given path
    MeritFun = zeros(size(im, 1), size(im, 2) + 1);        % allocate memory
    MeritFun(1, :) = ( 0:size(im, 2) ) * min(im(:)); % first row contains the minimum value that any path could add up to
    MeritFun(end, :) = MeritFun(1, :);% and also the final row contains this value
    costcont = zeros(size(im,1), size(im, 2));                               % allocate memory
    

    for j=1:size(im,2)             % over azimutal distance, 
            for m=2+vertical_freedom_pts:size(costcont, 1) - 1 - vertical_freedom_pts % over range distance/depth
                [A, n] = max(MeritFun(m - 1 - vertical_freedom_pts:m + 1 + vertical_freedom_pts, j));
                if n==1 && (MeritFun(m - 1, j)==MeritFun(m, j))
                    n=2;
                end
                MeritFun(m, j + 1) = im(m, j) + A;
                costcont(m, j) = m + (n - 2 - vertical_freedom_pts);
            end
    end
    
    
    % Trace back the optimal contour (that with the maximal cost) for periostal surface
    [~, n]=max(MeritFun(:, end)); % n is the depth-value where the cost function is lowest, so, this is the endpoint of our path
    Contour_Periosteum=zeros(1, size(im,2)); % allocate memory for the array of indices of depth, which form the path..
    Contour_Periosteum(end)=n;               % this is the end-point of the path...
    
    for j = (size(im,2)-1):-1:1 % Step through the azimutal distance to find the right path
            Contour_Periosteum(j) = costcont (Contour_Periosteum(j + 1), j + 1);
    end
    
    

% figure(5577)
% plot(Contour_Periosteum, 'm','linewidth',2');

max_z_index = Contour_Periosteum;
% title('select ROI for the fit of bone surface')
% RECT = getrect; %[xmin ymin width height]
% ROI_x_index = round(RECT(1)):round(RECT(1)+RECT(3));        
% max_x_index = ROI_x_index;
max_x_index = 1:size(im,2);
% max_z_index = max_z_index(max_x_index);