function [registeredPressure, V, LVP_average_unshifted] = main_HaemoFitting_LVP(directory, model_directory, dsMRI, edMRI, esMRI, eivcMRI, eivrMRI, outPressureFile)

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
%
% Outputs: 1) registeredPressure - registered pressure (DS = 0 kPa)
%          2) V - volumes over cardiac cycle
%          3) LVP_average_unishifted - registered pressure unshifted (DS != 0 kPa)
%
% Originally written by: Angus Cheng + Jenny Wang
%
% Adapted by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 26 June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Whether or not to recalculate haemodynamic points and choose cycles to use again
RECALC = 0;
chooseCycles = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save cardiac events from CIM

% Mat file to save cardiac events from MRI
cardiacEventsMRI_MatFile = sprintf('%s/cardiacEvents_MRI.mat', directory);
save(cardiacEventsMRI_MatFile, 'dsMRI', 'edMRI', 'esMRI', 'eivcMRI', 'eivrMRI');


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

% Plot raw LVP
FH = figure;
plot(LV_time, LV_Pressure, '-', 'markersize', 2);
title('Raw LV Pressure');
ylabel('Pressure (mmHg)');
xlabel('Time (ms)');

% Save figure
saveas(FH, sprintf('%s/raw-pressure-trace.png', directory));
close(FH)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert pressure to kPa

LV_Pressure = LV_Pressure / 7.50061561303;


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

subplot(2,1,1)
plot(LV_x_PS, 2*abs(LV_PS(1:LV_NFFT/2+1)));
title('LVP Power Spectrum')
ylabel('Amplitude')
xlabel('Frequency (Hz)')

% Plot zoomed in figure
subplot(2,1,2)
plot(LV_x_PS, 2*abs(LV_PS(1:LV_NFFT/2+1)));
title('LVP Power Spectrum - Zoomed')
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
    
    [cut_off_frequency, ~] = getpts(ftmp);
    close(ftmp);
    
    % Filter pressure
    LV_Pressure_filtered = filter_pressure(LV_Pressure, LV_samplerate(1), cut_off_frequency);
    
    % Plot filtered LVP
    FH3 = figure;
    plot(LV_time, LV_Pressure_filtered, '-', 'markersize', 2);
    title('Filtered LV Pressure');
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
        save(filterMatFile, 'cut_off_frequency', 'LV_Pressure');
    end
    
end

% If mat file exists, just load it with the filter cut off values
load(filterMatFile);

% Print to log file
fprintf('LV pressure trace low-pass cutoff frequency: %.3f\n', cut_off_frequency(1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find R peaks

% LVP
[LV_R_positions, LV_R_times, ~] = find_RR(LV_time, LV_ECG, LV_samplerate );

% Plot ECG traces with R peaks
FH4 = figure;
plot(LV_time, LV_ECG, '-', 'markersize', 2)
hold on
scatter(LV_R_times, LV_ECG(LV_R_positions), 5, 'ro', 'filled')
title('LV ECG - R Peaks');
ylabel('Voltage (mV)');
xlabel('Time (s)');

% Save figure
saveas(FH4, sprintf('%s/ecg-traces-r-peaks.png', directory));
close(FH4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose which LVP cycles to use

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
    
    save(cyclesMatFile, 'LV_cycles_to_keep');
    
else
    % If the mat file exists, load it
    load(cyclesMatFile);
    
end

% Print to log file
disp('LVP cycles to use:')
disp(LV_cycles_to_keep)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract all cardiac cycles

% Ouptut of extract_individual_cycles.m gives a cell array containing 
% time, pressure and ecg data for ONLY the cycles contained in LV_cycles_to_keep

% LVP
LV_cycles = extract_individual_cycles(LV_time, LV_R_positions, LV_cycles_to_keep, LV_Pressure, LV_ECG);

% Check: Overlay all LV pressure traces
FH5 = figure;
for i = 1:size(LV_cycles,2)
    plot(LV_cycles{1,i}, LV_cycles{2,i}, 'b-', 'markersize', 2);
    hold on;
end
title('All LV Pressure Traces')
xlabel('Time (s)')
ylabel('Pressure (kPa)')

% Save figure
saveas(FH5, sprintf('%s/all-cycles.png', directory));
close(FH5)

%%%%%%%%%%%%%%%%%%%%%%%% MARK CARDIAC EVENTS IN PRESSURE TRACES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

haemoEventMatFile = sprintf('%s/cardiacEvents_Haemo.mat', directory);

if ~exist(haemoEventMatFile, 'file') || RECALC
       
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
    save(haemoEventMatFile, 'ED', 'DS', 'eIVR', 'maxLVP', 'ED_t', 'DS_t', 'eIVR_t', 'maxLVP_t');
    
    
else
    load(haemoEventMatFile);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot cycles with identified landmarks

% Overlay all LV pressure traces
FH6 = figure;
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

% Save figure
saveas(FH6, sprintf('%s/all-cycles-landmarks.png', directory));
close(FH6)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get RR interval, MRI image timing and volumes

% Get mean RR interval from trigger times in image header files
%[RR_mean, TS_mean, TS_std, numFrames] = ExtractTriggerTime(model_directory);

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

timepoints_MRI = [edMRI, dsMRI, numFrames];

[LVP_average_unshifted, MRI_LVP_ed2ds, MRI_LVP_ds2ed] = InterpolatePressure_2Points(timepoints_MRI, LV_cycles, DS);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get final pressure trace to use as boundary conditions for parameter estimation

% Shift pressure so that LVP(DS) = 0
DS_pressure = min(LVP_average_unshifted); % Diastasis pressure
LVP_average = LVP_average_unshifted - DS_pressure; % Shifted LV pressure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create plots

plots_2points(timepoints_MRI, mri_time, LVP_average_unshifted, MRI_LVP_ed2ds, ...
        MRI_LVP_ds2ed, V, directory);


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
