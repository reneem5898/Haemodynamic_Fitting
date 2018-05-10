function [registeredPressure, V, LVP_average_unshifted] = main_HaemoFitting(directory, model_directory, dsMRI, edMRI, esMRI, eivcMRI, eivrMRI, outPressureFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function:
%
% Inputs: 1) directory - directory where pressure data are located and where
% all figures, output files will be saved (specific to volunteer/patient/animal)
%         2) model_directory - CIM model directory (specific to volunteer/
% patient/animal)
%         3) dsMRI - diastasis time point in MRI cine dataset
%         4) edMRI - end-diastolic time point in MRI cine dataset
%         5) esMRI - end-systolic time point in MRI cine dataset
%         6) eivcMRI - end isovolumetric contraction time point in MRI cine dataset
%         7) eivrMRI - end isovolumetric relaxation time point in MRI cine dataset
%         8) outPressureFile - filename for registered pressure file (THIS 
% FILE IS THE WHOLE POINT OF THIS CODE)

% Originally written by: Angus Cheng + Jenny Wang
%
% Adapted by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 20 February 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Whether or not to recalculate haemodynamic points and choose cycles to use again
RECALC = 0;
chooseCycles = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save cardiac events from CIM

% Mat file to save cardiac events from MRI
cardiacEventsMRI_MatFile = sprintf('%s/cardiacEvents_MRI.mat', directory);
save(cardiacEventsMRI_MatFile, 'dsMRI', 'edMRI', 'esMRI', 'eivcMRI', 'eivrMRI');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Get relevant cardiac events from MRI images - USER NEEDS TO ENTER THEM
% 
% % Mat file to save cardiac events from MRI
% cardiacEventsMRI_MatFile = sprintf('%s/cardiacEvents_MRI.mat', directory);
% 
% % If mat file doesn't exist
% if ~exist(cardiacEventsMRI_MatFile, 'file')
%     
%     % Get user to input relvant cardiac time points
%     prompt = {'End-diastole:', 'end IVC', 'End-systole', 'end IVR', 'Diastasis:'};
%     dlg_title = 'Cardiac Events: MRI';
%     num_lines = 1;
%     answer = inputdlg(prompt, dlg_title, num_lines);
%     
%     % Convert the string inputs to numbers
%     edMRI = str2double(answer{1});
%     eivcMRI = str2double(answer{2});
%     esMRI = str2double(answer{3});
%     eivrMRI = str2double(answer{4});
%     dsMRI = str2double(answer{5});
%     
%     % Save time points in mat file
%     save(cardiacEventsMRI_MatFile, 'dsMRI', 'edMRI', 'esMRI', 'eivcMRI', 'eivrMRI');
%     
% else
%     % Else if it exists, use previously chosen time points
%     load(cardiacEventsMRI_MatFile);
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract data - CHANGE THIS SECTION AS NEED TO READ IN YOUR DATA

% LVP and AOP file names
pressureFilenamesMat = sprintf('%s/pressureFileNames.mat', directory);
file_spec = sprintf('%s/*.mat', directory);
if ~exist(pressureFilenamesMat, 'file')
    LVP_file = uigetfile(file_spec, 'Choose LVP File');
    AOP_file = uigetfile(file_spec, 'Choose AOP File');
    save(pressureFilenamesMat, 'LVP_file', 'AOP_file');
end

load(pressureFilenamesMat);

% Load pressure/ecg data from Matlab file
% NOTE: for different input files, another function may be written in place
% of this load_data.m. OSU pressure/ECG data came in .mat files. However,
% for other studies, input may be in a different file type. Output should be:
% a) time vector (x-axis for pressure and ECG traces), b) pressure vector,
% c) ECG vector and d) samplerate (data points / second) [ 2 x 1] = [pressure_samplerate, ecg_samplerate]

% LVP
disp('Reading in LV pressure...');
disp(LVP_file);
[LV_time, LV_Pressure, LV_ECG, LV_samplerate] = load_data_LVP(directory, LVP_file);

% Hack - some pig data has different length ECG and pressure data - not sure what to do....
if length(LV_Pressure) ~= length(LV_ECG)
    
    disp('WARNING: LV Pressure and ECG traces are not the same length.')
    
    % Get vector of shortest length
    minLen = min([length(LV_Pressure), length(LV_ECG)]);
    
    % Update vectors so that all are the same length
    LV_ECG = LV_ECG(1:minLen);
    LV_Pressure = LV_Pressure(1:minLen);
    LV_time = LV_time(1:minLen);
end

% AOP
disp('Reading in AO pressure...');
disp(AOP_file);
[AO_time, AO_Pressure, AO_ECG, AO_samplerate] = load_data_AOP(directory, AOP_file);

% Hack - some pig data has different length ECG and pressure data - not sure what to do....
if length(AO_Pressure) ~= length(AO_ECG)
    
    disp('WARNING: Aortic Pressure and ECG traces are not the same length.')
    
    % Get vector of shortest length
    minLen = min([length(AO_Pressure), length(AO_ECG)]);
    
    % Update vectors so that all are the same length
    AO_ECG = AO_ECG(1:minLen);
    AO_Pressure = AO_Pressure(1:minLen);
    AO_time = AO_time(1:minLen);
end

% Plot raw LVP
FH = figure('Position', [100, 100, 600, 600]);
subplot(2,1,1)
plot(LV_time, LV_Pressure, '-', 'markersize', 2);
title('Raw LV Pressure');
ylabel('Pressure (mmHg)');
xlabel('Time (ms)');

% Plot raw AOP
subplot(2,1,2)
plot(AO_time, AO_Pressure, '-', 'markersize', 2);
title('Raw Aortic Pressure');
ylabel('Pressure (mmHg)');
xlabel('Time (ms)');

% Save figure
saveas(FH, sprintf('%s/raw-pressure-traces.png', directory));
close(FH)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert pressure to kPa

LV_Pressure = LV_Pressure / 7.50061561303;
AO_Pressure = AO_Pressure / 7.50061561303;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create power spectrum - pressure signals
% Plot to allow user to choose cut-off frequencies for low-pass filters

% Plot power spectrum
FH2 = figure('Position', [100, 100, 1200, 600]);

%%%%%%%%%%%%%%%%%%%%%% LVP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform - LVP
L = length(LV_Pressure);
LV_NFFT = 2^nextpow2(L);
LV_PS = fft(LV_Pressure, LV_NFFT)/L;
LV_x_PS = LV_samplerate(1)/2*linspace(0,1,LV_NFFT/2+1);

subplot(2,2,1)
plot(LV_x_PS, 2*abs(LV_PS(1:LV_NFFT/2+1)));
title('LVP Power Spectrum')
ylabel('Amplitude')
xlabel('Frequency (Hz)')

% Plot zoomed in figure
subplot(2,2,3)
plot(LV_x_PS, 2*abs(LV_PS(1:LV_NFFT/2+1)));
title('LVP Power Spectrum - Zoomed')
ylabel('Amplitude')
xlabel('Frequency (Hz)')
xlim([0 25])

%%%%%%%%%%%%%%%%%%%%%% AOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier transform - AOP
L = length(AO_Pressure);
AO_NFFT = 2^nextpow2(L);
AO_PS = fft(AO_Pressure, AO_NFFT)/L;
AO_x_PS = AO_samplerate(1)/2*linspace(0,1,AO_NFFT/2+1);

subplot(2,2,2)
plot(AO_x_PS, 2*abs(AO_PS(1:AO_NFFT/2+1)));
title('AOP Power Spectrum')
ylabel('Amplitude')
xlabel('Frequency (Hz)')

% Plot zoomed in figure
subplot(2,2,4)
plot(AO_x_PS, 2*abs(AO_PS(1:AO_NFFT/2+1)));
title('AOP Power Spectrum - Zoomed')
ylabel('Amplitude')
xlabel('Frequency (Hz)')
xlim([0 25])

% Save figure
saveas(FH2, sprintf('%s/pressure-power-spectra.png', directory));
close(FH2)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Low-pass filtering

% Save low-pass filter cut off values in a mat file (to avoid repeating this step if running again)
filterMatFile = sprintf('%s/filterValues.mat', directory);

while ~exist(filterMatFile, 'file')
    
    %%%%%%%%%%%%%%%%%%% LVP %%%%%%%%%%%%%%%%%%%%%%%%%%
    waitfor(msgbox('Select the cutoff frequency for low-pass filtering the LV pressure trace. 1) Click on the graph at the frequency (x location) to put cut-off, then 2) Press enter.', 'Cut-off Frequency'));
    
    % Get user to select point on the graph for low-pass filter cut off frequency
    ftmp = figure;
    plot(LV_x_PS, 2*abs(LV_PS(1:LV_NFFT/2+1)));
    title('LVP Power Spectrum - Zoomed', 'FontSize', 16)
    ylabel('Amplitude')
    xlabel('Frequency (Hz)')
    xlim([0 25])
    
    [cut_off_frequency(1), ~] = getpts(ftmp);
    close(ftmp);
    
    %%%%%%%%%%%%%%%%%%% AOP %%%%%%%%%%%%%%%%%%%%%%%%%%
    waitfor(msgbox('Select the cutoff frequency for low-pass filtering the AORTIC pressure trace. 1) Click on the graph at the frequency (x location) to put cut-off, then 2) Press enter.', 'Cut-off Frequency'));
    
    % Get user to select point on the graph for low-pass filter cut off frequency
    ftmp = figure;
    plot(AO_x_PS, 2*abs(AO_PS(1:AO_NFFT/2+1)));
    title('AOP Power Spectrum - Zoomed', 'FontSize', 16)
    ylabel('Amplitude')
    xlabel('Frequency (Hz)')
    xlim([0 25])
    
    [cut_off_frequency(2), ~] = getpts(ftmp);
    close(ftmp);
    
    % % Get cut off frequencies for filtering
    % prompt = {'Enter cut-off frequency for low-pass filter: LVP', 'Enter cut-off frequency for low-pass filter: AOP'};
    % dlg_title = 'Input';
    % num_lines = 1;
    % defaultans = {'25', '25'};
    % cut_off_frequency = inputdlg(prompt, dlg_title, num_lines, defaultans);
    
    % LV
    %[LV_Pressure] = filter_pressure(LV_Pressure, LV_samplerate(1), str2double(cut_off_frequency{1}));
    LV_Pressure_filtered = filter_pressure(LV_Pressure, LV_samplerate(1), cut_off_frequency(1));
    %LV_Pressure = medfilt1(LV_Pressure, 20);
    
    % AOP
    %[AO_Pressure] = filter_pressure(AO_Pressure, AO_samplerate(1), str2double(cut_off_frequency{2}));
    AO_Pressure_filtered = filter_pressure(AO_Pressure, AO_samplerate(1), cut_off_frequency(2));
    %AO_Pressure = medfilt1(AO_Pressure, 20);
    
    % Plot filtered LVP
    FH3 = figure('Position', [100, 100, 600, 600]);
    subplot(2,1,1)
    plot(LV_time, LV_Pressure_filtered, '-', 'markersize', 2);
    title('Filtered LV Pressure');
    ylabel('Pressure (kPa)');
    xlabel('Time (s)');
    
    % Plot filtered AOP
    subplot(2,1,2)
    plot(AO_time, AO_Pressure_filtered, '-', 'markersize', 2);
    title('Filtered Aortic Pressure_filtered');
    ylabel('Pressure (kPa)');
    xlabel('Time (s)');
    
    % Save figure
    saveas(FH3, sprintf('%s/filtered-pressure-traces.png', directory));
    
    % Ask user if they are satisfied with the smoothing - if no, repeat process of choosing cut-off frequency
    qstring = 'Are you satisfied with the smoothing?';
    tmp = questdlg(qstring,'Finished Smoothing?','Yes','No','Yes');
    if strfind(tmp, 'Yes')
        
        % Save smoothed pressure traces
        LV_Pressure = LV_Pressure_filtered;
        AO_Pressure = AO_Pressure_filtered;
        save(filterMatFile, 'cut_off_frequency', 'AO_Pressure', 'LV_Pressure');
    end
    
end

% If mat file exists, just load it with the filter cut off values
load(filterMatFile);


% Print to log file
fprintf('LV pressure trace low-pass cutoff frequency: %.3f\n', cut_off_frequency(1));
fprintf('Aortic pressure trace low-pass cutoff frequency: %.3f\n', cut_off_frequency(2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find R peaks

% LVP
[LV_R_positions, LV_R_times, LV_RR] = find_RR(LV_time, LV_ECG, LV_samplerate );
% AOP
[AO_R_positions, AO_R_times, AO_RR] = find_RR(AO_time, AO_ECG, AO_samplerate );

% Plot ECG traces with R peaks
FH4 = figure('Position', [100, 100, 600, 600]);
subplot(2,1,1)
plot(LV_time, LV_ECG, '-', 'markersize', 2)
hold on
scatter(LV_R_times, LV_ECG(LV_R_positions), 5, 'ro', 'filled')
title('LV ECG - R Peaks');
ylabel('Voltage (mV)');
xlabel('Time (s)');

subplot(2,1,2)
plot(AO_time, AO_ECG, '-', 'markersize', 2)
hold on
scatter(AO_R_times, AO_ECG(AO_R_positions), 5, 'ro', 'filled')
title('Aortic ECG - R Peaks');
ylabel('Voltage (mV)');
xlabel('Time (s)');

% Save figure
saveas(FH4, sprintf('%s/ecg-traces-r-peaks.png', directory));
close(FH4)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose which LVP and AOP cycles to use

% Description of how to choose cycles to use: Jenny?
% Jenny - in St. Francis dataset, choose cycles that look like there is a pattern. There is often
% a pattern due to respiratory wave. In this case, choose the higher pressure cycles. It's okay
% if only a few cycles are chosen.
% HOWEVER, for the pig data, I am choosing most all cycles and only leaving 
% out ones which are clearly different. For example, in some datasets, the ECG
% doesn't seem to be gated and there are two pressure cycles for one R peak. That cycle
% has been excluded since it is clearly not like the rest. In the pig data, there
% are numerous cycles and so any differences between individual cycles will
% be averaged out. 
% In general, how to choose which cycles to use and which to leave out will
% be at the discretion of the user. 

% Mat file
cyclesMatFile = sprintf('%s/cyclesToUse.mat', directory);

% If mat file doesn't exist, get user to select which cycles to use for analysis
if ~exist(cyclesMatFile, 'file') || chooseCycles
    
    % Choose LVP cycles to keep
    LV_cycles_fig = sprintf('%s/LVP_cycles.png', directory);
    LV_cycles_to_keep = keepCycles(LV_Pressure, LV_R_positions, LV_cycles_fig); % LV trace
    
    % Choose AOP cycles to keep
    AO_cycles_fig = sprintf('%s/AOP_cycles.png', directory);
    AO_cycles_to_keep = keepCycles(AO_Pressure, AO_R_positions, AO_cycles_fig); % Aortic trace
    
    save(cyclesMatFile, 'LV_cycles_to_keep', 'AO_cycles_to_keep');
    
else
    % If the mat file exists, load it
    load(cyclesMatFile);
    
end

% Print to log file
disp('LVP cycles to use:')
disp(LV_cycles_to_keep)
disp('AOP cycles to use:')
disp(AO_cycles_to_keep)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract all cardiac cycles

% Ouptut of extract_individual_cycles.m gives a cell array containing 
% time, pressure and ecg data for ONLY the cycles contained in LV_cycles_to_keep
% and AO_cycles_to_keep

% LVP
LV_cycles = extract_individual_cycles(LV_time, LV_R_positions, LV_cycles_to_keep, LV_Pressure, LV_ECG);
% AOP
AO_cycles = extract_individual_cycles(AO_time, AO_R_positions, AO_cycles_to_keep, AO_Pressure, AO_ECG);

% Check: Overlay all LV pressure traces
FH5 = figure('Position', [100, 100, 600, 600]);
subplot(2,1,1)
for i = 1:size(LV_cycles,2)
    plot(LV_cycles{1,i}, LV_cycles{2,i}, 'b-', 'markersize', 2);
    hold on;
end
title('All LV Pressure Traces')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
hold off;

% Check: Overlay all AO pressure traces
subplot(2,1,2)
for i = 1:size(AO_cycles,2)
    plot(AO_cycles{1,i}, AO_cycles{2,i}, 'b-', 'markersize', 2);
    hold on;
end
title('All Aortic Pressure Traces')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
hold off;

% Save figure
saveas(FH5, sprintf('%s/all-cycles.png', directory));
close(FH5)

%%%%%%%%%%%%%%%%%%%%%%%% MARK CARDIAC EVENTS IN PRESSURE TRACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

haemoEventMatFile = sprintf('%s/cardiacEvents_Haemo.mat', directory);

if ~exist(haemoEventMatFile, 'file') || RECALC
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get eIVC and ES from AOP traces
    
    % Initialise variables
    eIVC = zeros(size(AO_cycles,2),1);
    eIVC_t = zeros(size(AO_cycles,2),1);
    ES = zeros(size(AO_cycles,2),1);
    ES_t = zeros(size(AO_cycles,2),1);
    
    % Directions for the user
    waitfor(msgbox('Select the peak which corresponds to the dicrotic notch. If it was not identified, choose NaN.', 'Dicrotic notch'));
    
    % Loop through each cycle in the aortic trace
    for i = 1:size(AO_cycles,2)
        
        % End IVC = min(AOP)
        [~,l] = findpeaks(-AO_cycles{2,i}); % Find local minima in AOP
        eIVC(i) = l(1); % Get first local minimum
        %eIVC(i) = find(AO_cycles{2,i} == min(AO_cycles{2,i}));
        
        % ES = dicrotic notch in aortic pressure trace
        %%%%% REQUIRES USER INPUT %%%%%
        ES(i) = getES_dicroticNotch(AO_cycles{2,i});
        
        % Get times corresponding to eIVC and ES
        eIVC_t(i) = AO_cycles{1,i}(eIVC(i));
        ES_t(i) = AO_cycles{1,i}(ES(i));
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get eIVR, diastasis and ED from LVP traces
    
    % Initialise variables
    DS = zeros(size(LV_cycles,2),1);
    ED = zeros(size(LV_cycles,2),1);
    eIVR = zeros(size(LV_cycles,2),1);
    maxLVP = zeros(size(LV_cycles,2),1);
    DS_t = zeros(size(LV_cycles,2),1);
    ED_t = zeros(size(LV_cycles,2),1);
    eIVR_t = zeros(size(LV_cycles,2),1);
    maxLVP_t = zeros(size(LV_cycles,2),1);
    
    % Loop through each cycle in the LV trace
    for i = 1:size(LV_cycles,2)
        
        % ED = 1
        ED(i) = 1; % ED (R-peak) is always first point in each cycle (see extract_individual_cycles.m)
        
        % Calculate second derivative of LVP
        dLVP_dt = gradient(LV_cycles{2,i}, 1/LV_samplerate(1));
        d2LVP_dt2 = gradient(dLVP_dt, 1/LV_samplerate(1));
        
        % Diastasis
        middleIndex = floor(length(LV_cycles{2,i})/2);
        DS(i) = find(LV_cycles{2,i} == min(LV_cycles{2,i}(middleIndex:end))); % Minimum pressure in second half of LVP trace
        %DS(i) = getDS(LV_cycles{2,i});
        
        % eIVR = POI in LVP -- RM: need citation for this.
        eIVR(i) = get_eIVR(LV_cycles{2,i}, d2LVP_dt2);
        
        % Peak LV pressure
        maxLVP(i) = find(LV_cycles{2,i} == max(LV_cycles{2,i}));
        
        % Get times corresponding to ED, eIVR and DS
        ED_t(i) = LV_cycles{1,i}(ED(i));
        DS_t(i) = LV_cycles{1,i}(DS(i));
        eIVR_t(i) = LV_cycles{1,i}(eIVR(i));
        maxLVP_t(i) = LV_cycles{1,i}(maxLVP(i));
        
    end
    
    % Save cardiac events in .mat file
    save(haemoEventMatFile, 'ES', 'ED', 'DS', 'eIVR', 'eIVC', 'maxLVP', 'ES_t', 'ED_t', 'DS_t', 'eIVR_t', 'eIVC_t', 'maxLVP_t');
    
    
else
    load(haemoEventMatFile);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot cycles with identified landmarks

% Overlay all LV pressure traces
FH6 = figure('Position', [100, 100, 600, 600]);
subplot(2,1,1)
for i = 1:size(LV_cycles,2)
    plot(LV_cycles{1,i}, LV_cycles{2,i}, 'b-', 'markersize', 2);
    hold on;
    scatter(LV_cycles{1,i}(DS(i)), LV_cycles{2,i}(DS(i)), 'ro', 'filled')
    scatter(LV_cycles{1,i}(ED(i)), LV_cycles{2,i}(ED(i)), 'ro', 'filled')
    scatter(LV_cycles{1,i}(eIVR(i)), LV_cycles{2,i}(eIVR(i)), 'ro', 'filled')
end
title('All LV Pressure Traces with Landmarks')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
hold off;

% Overlay all AO pressure traces
subplot(2,1,2)
for i = 1:size(AO_cycles,2)
    plot(AO_cycles{1,i}, AO_cycles{2,i}, 'b-', 'markersize', 2);
    hold on;
    scatter(AO_cycles{1,i}(eIVC(i)), AO_cycles{2,i}(eIVC(i)), 'ro', 'filled')
    scatter(AO_cycles{1,i}(ES(i)), AO_cycles{2,i}(ES(i)), 'ro', 'filled')
end
title('All Aortic Pressure Traces with Landmarks')
xlabel('Time (s)')
ylabel('Pressure (kPa)')
hold off;

% Save figure
saveas(FH6, sprintf('%s/all-cycles-landmarks.png', directory));
close(FH6)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get RR interval, MRI image timing and volumes

% Get mean RR interval from trigger times in image header files
%[RR_mean, TS_mean, TS_std, numFrames] = ExtractTriggerTime(MR_directory);

% Get RR interval from CIM model file - IF NO TRIGGER TIMES ARE AVAILABLE IN HEADER, USE THIS LINE INSTEAD :-)
[RR_mean, numFrames] = getRR_CIM(model_directory);
mri_time = linspace(1,RR_mean,numFrames)/1000;

% Get LV model volumes
V = getModelVolumes(model_directory, numFrames);

% If ED is not the first frame, make it the first frame (fix volumes)
if edMRI ~= 1

    % Save old time points and old volumes for shifting data back
    timePointsOld = [edMRI, eivcMRI, esMRI, eivrMRI, dsMRI];
    V_Old = V;
    
    % Shifted volumes
    V = [V(edMRI:end); V(1:(edMRI-1))];
    
    d = edMRI - 1; % Get difference
    edMRI = 1; % Reset ed to 1
    
    % Shift all other cardiac time points
    eivcMRI = eivcMRI - d;
    esMRI = esMRI - d;
    eivrMRI = eivrMRI - d;
    dsMRI = dsMRI - d;
end
    

% % Check to see if initial and final volumes are within 2% of each other (2% chosen arbitrarily - can change)
% V_error = abs(V(end) - V(1)); % Error between first and last frames
% 
% if V_error > 0.02 * V(1) % If error is greater than 2%
%    V = [V(end) V]; % Replicate last frame
%    numFrames = numFrames + 1;
%    
%    disp('WARNING: LV volumes were significantly different between last and first frame. The last frame was added to the beginning of the cycle.')
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate pressures at MRI frames 

num_MRI_frames_ed2eivc = eivcMRI - edMRI + 1; % Number of frames in images from ED to end IVC 
num_MRI_frames_eivc2es = esMRI - eivcMRI + 1; % Number of frames in images from end IVC to ES 
num_MRI_frames_es2eivr = eivrMRI - esMRI + 1; % Number of frames in images from ES to end IVR 
num_MRI_frames_eivr2ds = dsMRI - eivrMRI + 1; % Number of frames in images from end IVR to DS 
num_MRI_frames_ds2ed = numFrames - dsMRI + edMRI; % Number of frames in images from DS to ED (excluding last point)

c = 1; % counter
eIVC_LV = zeros(size(LV_cycles,2)*size(AO_cycles,2),1); % Initialise variable to save all possible eIVC points
ES_LV = zeros(size(LV_cycles,2)*size(AO_cycles,2),1); % Initialise variable to save all possible ES indices
%ESt_LV = zeros(size(LV_cycles,2)*size(AO_cycles,2),1); % Initialise variable to save all possible ES time points

% Initialise variables to save pressures at MRI frames
% Sizes differ depending on options
MRI_LVP_ed2eivc = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_ed2eivc); 
MRI_LVP_eivc2es = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_eivc2es-1); % -1 because data will be exclusive of eivc
MRI_LVP_es2eivr = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_es2eivr-1); % -1 because data will be exclusive of es
MRI_LVP_eivr2ds = zeros(size(LV_cycles,2), num_MRI_frames_eivr2ds-2); % -1 because data will be exclusive of eivr AND ds
MRI_LVP_ds2ed = zeros(size(LV_cycles,2), num_MRI_frames_ds2ed); 

for i = 1:size(LV_cycles,2)
    
    LVP_c = LV_cycles{2,i}; % Current LV pressure cycle
    
    LVP_ed2MaxLVP = LVP_c(ED(i):maxLVP(i)); % Section from ED to maximum LVP
    LVP_maxLVP2ds = LVP_c(maxLVP(i):DS(i)); % Section from maximum LVP to DS
    
    %% eIVR to DS
    
    % Interpolate pressure at each MRI frame between ds and ed
    pressure_eivr2ds = LVP_c(eIVR(i):DS(i)); % LV pressure between DS and ED
    interpPressure = LinearInterpolatePressure(pressure_eivr2ds, length(pressure_eivr2ds), num_MRI_frames_eivr2ds);
    MRI_LVP_eivr2ds(i,:) = interpPressure(2:end-1); % Exclude end IVR
    
    %% DS to ED
    
    % Interpolate pressure at each MRI frame between ds and ed
    pressure_ds2ed = LVP_c(DS(i):end); % LV pressure between DS and ED
    interpPressure = LinearInterpolatePressure(pressure_ds2ed, length(pressure_ds2ed), num_MRI_frames_ds2ed);
    MRI_LVP_ds2ed(i,:) = interpPressure;
    
    for j = 1:size(AO_cycles,2)
        
        AOP_c = AO_cycles{2,j}; % Current AO pressure cycle
        
        %% ED to eIVC
        
        % Get point in LVP when LVP == AOP(eIVC)
        AOP_eIVC = AOP_c(eIVC(j)); % Aortic pressure at endIVC
        eIVC_LV(c) = find(abs(LVP_ed2MaxLVP - AOP_eIVC) == min(abs(LVP_ed2MaxLVP - AOP_eIVC)));
        
        % Interpolate pressure at MRI image frames from ED to eIVC
        pressure_ed2eivc = LVP_c(ED(i):eIVC_LV(c));
        interpPressure = LinearInterpolatePressure(pressure_ed2eivc, length(pressure_ed2eivc), num_MRI_frames_ed2eivc);
        MRI_LVP_ed2eivc(c,:) = interpPressure;
        
        %% eIVC to ES
        
        % Get point in LVP when LVP == AOP(ES)
        AOP_ES = AOP_c(ES(j)); % Aortic pressure at ES
        ES_LV(c) = find(abs(LVP_maxLVP2ds - AOP_ES) == min(abs(LVP_maxLVP2ds - AOP_ES))) + maxLVP(i);
        %ESt_LV(c) = LV_cycles{1,i}(ES_LV(c)) * LV_cycles{1,i}(DS(i))/mean(LV_cycles{1,i}(DS(:))); % DON'T UNDERSTAND THIS LINE - WHY SCALE BY DS TIMING????
        
        % Interpolate pressure at MRI image frames from eIVC to ES
        pressure_eivc2es = LVP_c(eIVC_LV(c):ES_LV(c));
        interpPressure = LinearInterpolatePressure(pressure_eivc2es, length(pressure_eivc2es), num_MRI_frames_eivc2es);
        MRI_LVP_eivc2es(c,:) = interpPressure(2:end); % Exclude eIVC
        
        %% ES to eIVR
        
        pressure_es2eivr = LVP_c(ES_LV(c):eIVR(i));
        interpPressure = LinearInterpolatePressure(pressure_es2eivr, length(pressure_es2eivr), num_MRI_frames_es2eivr);
        MRI_LVP_es2eivr(c,:) = interpPressure(2:end); % Exclude ES
        
        % Count up
        c = c + 1;
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get final pressure trace to use as boundary conditions for parameter estimation

% Average pressures and put together segments
%                   ED to eIVC                 eIVC to ES             ES to eIVR                 eIVR to DS               DS to ED
LVP_average_unshifted = [mean(MRI_LVP_ed2eivc,1), mean(MRI_LVP_eivc2es,1), mean(MRI_LVP_es2eivr,1), mean(MRI_LVP_eivr2ds,1), mean(MRI_LVP_ds2ed,1)]';

% Shift pressure so that LVP(DS) = 0
DS_pressure = min(LVP_average_unshifted); % Diastasis pressure
LVP_average = LVP_average_unshifted - DS_pressure; % Shifted LV pressure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create plots

%%%%%%%%%%%%%%%% Beat-averaged left ventricular pressure %%%%%%%%%%%%%%%%%%
FH8 = figure;
% Plot ED to eIVC
t = repmat(mri_time(edMRI:eivcMRI),size(MRI_LVP_ed2eivc,1),1);
scatter(t(:), MRI_LVP_ed2eivc(:), 'g*')
hold on
% Plot eIVC to ES
t = repmat(mri_time((eivcMRI+1):esMRI),size(MRI_LVP_eivc2es,1),1);
scatter(t(:), MRI_LVP_eivc2es(:), 'r*')
% Plot ES to eIVR
t = repmat(mri_time((esMRI+1):eivrMRI),size(MRI_LVP_es2eivr,1),1);
scatter(t(:), MRI_LVP_es2eivr(:), 'm*')
% Plot eIVR to DS
t = repmat(mri_time((eivrMRI+1):dsMRI-1),size(MRI_LVP_eivr2ds,1),1);
scatter(t(:), MRI_LVP_eivr2ds(:), 'c*')
% Plot DS to ED
if edMRI == 1
    t = repmat(mri_time(dsMRI:end),size(MRI_LVP_ds2ed,1),1);
else
    t = repmat([mri_time(dsMRI:end),mri_time(1:(edMRI-1))],size(MRI_LVP_ds2ed,1),1);
end
scatter(t(:), MRI_LVP_ds2ed(:), 'b*')
% Plot average
plot(mri_time, LVP_average_unshifted, 'k*-')
xlabel('Time (s)', 'FontSize', 16)
ylabel('LV Pressure (kPa)', 'FontSize', 16)
title('Beat-averaged LV Pressure', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH8, sprintf('%s/beat-averaged-lv-pressure.png', directory));


%%%%%%%%%%%%%%%%%%% Pressure-volume loop %%%%%%%%%%%%%%%%%%
FH9 = figure;
% Plot ED to eIVC
t = repmat(V(edMRI:eivcMRI)',size(MRI_LVP_ed2eivc,1),1);
scatter(t(:), MRI_LVP_ed2eivc(:), 'g*')
hold on
% Plot eIVC to ES
t = repmat(V((eivcMRI+1):esMRI)',size(MRI_LVP_eivc2es,1),1);
scatter(t(:), MRI_LVP_eivc2es(:), 'r*')
% Plot ES to eIVR
t = repmat(V((esMRI+1):eivrMRI)',size(MRI_LVP_es2eivr,1),1);
scatter(t(:), MRI_LVP_es2eivr(:), 'm*')
% Plot eIVR to DS
t = repmat(V((eivrMRI+1):dsMRI-1)',size(MRI_LVP_eivr2ds,1),1);
scatter(t(:), MRI_LVP_eivr2ds(:), 'c*')
% Plot DS to ED
if edMRI == 1
    t = repmat(V(dsMRI:end)',size(MRI_LVP_ds2ed,1),1);
else
    t = repmat([V(dsMRI:end),V(1:(edMRI-1))]',size(MRI_LVP_ds2ed,1),1);
end
scatter(t(:), MRI_LVP_ds2ed(:), 'b*')
plot(V, LVP_average_unshifted, 'k*-')
xlabel('Volume (mL)', 'FontSize', 16)
ylabel('LV Pressure (kPa)', 'FontSize', 16)
title('Pressure-Volume Loop', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH9, sprintf('%s/pv-loop.png', directory));


%%%%%%%%%%%%%%%%%%% Volume with cardiac events %%%%%%%%%%%%%%%%%%
FH10 = figure;
plot(mri_time, V, 'ksq-', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold on
scatter(mri_time(dsMRI), V(dsMRI), 50, 'ro', 'filled')
scatter(mri_time(edMRI), V(edMRI), 50, 'ro', 'filled')
scatter(mri_time(esMRI), V(esMRI), 50, 'ro', 'filled')
scatter(mri_time(eivrMRI), V(eivrMRI), 50, 'ro', 'filled')
scatter(mri_time(eivcMRI), V(eivcMRI), 50, 'ro', 'filled')
xlabel('Time (s)', 'FontSize', 16)
ylabel('Volume (mL)', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH10, sprintf('%s/volume-mri-cardiac-events.png', directory));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check to see if MRI points were shifted. If so, shift pressure to be
% aligned with MRI time points

if exist('timePointsOld', 'var')
    shiftTime = timePointsOld(1) - edMRI;
    LVP_average = [LVP_average((numFrames - shiftTime + 1):end); LVP_average(1:(numFrames-shiftTime))]; % Shift average pressures to be aligned with MRI time points
    V = V_Old; % Use old volumes
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write out final LVP at MRI time points to text file

% Create matrix: [ image frame number  LVP ]
registeredPressure = [linspace(1,length(V),length(V))' LVP_average];
dlmwrite(outPressureFile, registeredPressure, 'delimiter', '\t', 'precision', '%.5g');
