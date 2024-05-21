function [acquisition, reconstruction] = GenerateParamRecon(recorded, parameters)
    C_LENS = 1540; LENS_THICKNESS = 1e-6;
    
    RAW_RF_SIG = recorded.Signals;
    acquisition.Nsamples = recorded.NbTimePoints; % Number of time samples
    acquisition.Fs = 1/recorded.Temporal_step_us*1e6; % in Hz
    acquisition.PITCH = recorded.Pitch*recorded.Spatial_step_mm*1e-3;
    acquisition.NELEMENTS = recorded.NbElements; % Number of receivers
    acquisition.NSOURCES = size(RAW_RF_SIG,3);   % Number of emissions
    acquisition.FREQ_Transducteur = 2.5e6;
    acquisition.ELEM_WIDTH = recorded.Width*recorded.Spatial_step_mm*1e-3;
    
    acquisition.XS = (0:acquisition.NSOURCES-1)*acquisition.PITCH;
    acquisition.XS = acquisition.XS - mean(acquisition.XS);
    
    
    %
    pixelSize = 1/4; % 1/4 de longueur d'onde
    
    xmin = -10e-3;%
    xmax = -xmin;%
    zmin = 2e-3;%
    zmax = 15e-3;%
    
    acquisition.ZS = zeros(size(acquisition.XS));
    acquisition.C_LENS = C_LENS;
    acquisition.LENS_THICKNESS = LENS_THICKNESS;
    field_probe_struct = {'Fs','FREQ_Transducteur','PITCH',...
    'C_LENS','LENS_THICKNESS','NELEMENTS','Nsamples',...
    'XS','ZS','NSOURCES','offset'};
    nb_wavelengths = 3;
    input_signal_duration = nb_wavelengths/acquisition.FREQ_Transducteur;
    acquisition.offset = round(input_signal_duration/2*acquisition.Fs);
    assert(all(isfield(acquisition,field_probe_struct)))
    
    %
    BEST_C_LAYER_2 = 3500;
    if ~isfield(parameters.interface, 'periost') || parameters.interface.periost == 0
        BEST_C_LAYER_1 = BEST_C_LAYER_2;
    else 
        BEST_C_LAYER_1 = 1540;
    end  
    %
    D_PROBE_PERIOS = 10e-3;
    tx_fnumber = 0.6; % Transmit fnumber
    rx_fnumber = 0.6; %  Receive fnumber
    reconstruction = struct('MIN_CORTICAL',3e-3, 'MAX_CORTICAL',15e-3,...
        'D_PROBE_PERIOS',D_PROBE_PERIOS,...
        'C_TISSUE',BEST_C_LAYER_1,'C_RADIAL',BEST_C_LAYER_2,'C_AXIAL',BEST_C_LAYER_2,...
        'C_MARROW',1540,'ANISO_SHAPE_COEF',0,...
        'FnumberMin',rx_fnumber, 'SubApertureApodis',2, ...
        'NeedTravelTime',0,'ReConTo',2,...
        'Max_err_allowed_ReceiveAngle_rad', deg2rad(5));
    TxHalfOpeningAngInSkinDeg=atand(1/2/tx_fnumber);
    reconstruction.HalfOpeningAngInLensRad = asin(acquisition.C_LENS/reconstruction.C_TISSUE*sind(TxHalfOpeningAngInSkinDeg)); % [rad]

    reconstruction.pixel_size = pixelSize*BEST_C_LAYER_1/acquisition.FREQ_Transducteur;
    reconstruction.X = xmin:reconstruction.pixel_size:xmax;
    reconstruction.Z = zmin:reconstruction.pixel_size:zmax;
    reconstruction.Xmm = reconstruction.X*1e3; reconstruction.Zmm = reconstruction.Z*1e3;

    % assert(all(isfield(reconstruction,field_recon_struc)))
    NCORE = 24;
    reconstruction.rf_data = RAW_RF_SIG;
    reconstruction.PROBE_PARAM = acquisition;
    reconstruction.Setup = [reconstruction.ANISO_SHAPE_COEF LENS_THICKNESS acquisition.PITCH...
    reconstruction.D_PROBE_PERIOS reconstruction.MIN_CORTICAL reconstruction.MAX_CORTICAL ...
    reconstruction.HalfOpeningAngInLensRad rx_fnumber...
    reconstruction.Max_err_allowed_ReceiveAngle_rad acquisition.Fs ...
    acquisition.C_LENS reconstruction.C_TISSUE reconstruction.C_MARROW...
    reconstruction.C_AXIAL reconstruction.C_RADIAL NCORE ...
    reconstruction.ReConTo reconstruction.NeedTravelTime reconstruction.SubApertureApodis];
    
    
    acquisition.XR = acquisition.XS;
    acquisition.ZR = acquisition.ZS;

    % TIME OF FLIGHT and ANGLES OF VIEW in the two cases  
    if  ~isfield(parameters.interface, 'periost') || parameters.interface.periost == 0
        [Angle_T,Angle_R] = function_get_angle_of_view_sa(reconstruction.X,reconstruction.Z,acquisition.XS,acquisition.ZS,...
        acquisition.XR,acquisition.ZR);
        [Time_T, Time_R] = function_get_time_of_flight_sa(reconstruction.X,reconstruction.Z,BEST_C_LAYER_1, ...
            acquisition.XS,acquisition.ZS, acquisition.XR,acquisition.ZR);
        
        reconstruction.Bone.timeFlight.Time_T =Time_T;
        reconstruction.Bone.timeFlight.Time_R = Time_R;
        reconstruction.Bone.Degre.angleView.Angle_T = Angle_T;
        reconstruction.Bone.Degre.angleView.Angle_R = Angle_R;
        reconstruction.Bone.Pixel.angleView.Angle_T = Angle_T;
        reconstruction.Bone.Pixel.angleView.Angle_R = Angle_R;
    
    else 
        [~, reconstruction.PERI, reconstruction.ENDO, Angles,Times] = DAS_boneimage_Salome(reconstruction, reconstruction.FnumberMin, ...
        xmin,xmax, false);
        % Time of flight
        reconstruction.Tissu.timeFlight.Time_T = Times.Time_T_Tissue;
        reconstruction.Tissu.timeFlight.Time_R = Times.Time_R_Tissue;

        reconstruction.Bone.timeFlight.Time_T = Times.Time_T_Bone;
        reconstruction.Bone.timeFlight.Time_R = Times.Time_R_Bone;
        
        % Angle of view in pixel
        reconstruction.Tissu.Pixel.angleView.Angle_T = rad2deg(Angles.APix_T_Tissue);
        reconstruction.Tissu.Pixel.angleView.Angle_R = rad2deg(Angles.APix_T_Tissue);
        reconstruction.Bone.Pixel.angleView.Angle_T = rad2deg(Angles.APix_T_Bone);
        reconstruction.Bone.Pixel.angleView.Angle_R = rad2deg(Angles.APix_T_Bone);

        % Angle of view in degree
        reconstruction.Bone.Degre.angleView.Angle_T = rad2deg(Angles.Angle_T_Bone);
        reconstruction.Bone.Degre.angleView.Angle_R = rad2deg(Angles.Angle_T_Bone);
        reconstruction.Tissu.Degre.angleView.Angle_T = rad2deg(Angles.Angle_T_Tissue);
        reconstruction.Tissu.Degre.angleView.Angle_R = rad2deg(Angles.Angle_R_Tissue);
    end 
end
