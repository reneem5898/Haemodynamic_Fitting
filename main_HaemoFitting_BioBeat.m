function [registeredPressure, V, LVP_average_unshifted] = main_HaemoFitting_BioBeat(directory, model_directory, dsMRI, edMRI, esMRI, eivcMRI, eivrMRI, outPressureFile)

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

% Hard coded for BioBeat
SAMPLE_RATE = 240;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Save cardiac events from CIM
% 
% % Mat file to save cardiac events from MRI
% cardiacEventsMRI_MatFile = sprintf('%s/cardiacEvents_MRI.mat', directory);
% save(cardiacEventsMRI_MatFile, 'dsMRI', 'edMRI', 'esMRI', 'eivcMRI', 'eivrMRI');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract data - CHANGE THIS SECTION AS NEED TO READ IN YOUR DATA

% LVP and AOP file names
pressureFilenamesMat = sprintf('%s/pressureFileNames.mat', directory);
file_spec = sprintf('%s/*.txt', directory);
if ~exist(pressureFilenamesMat, 'file')
    LVP_file = uigetfile(file_spec, 'Choose LVP File'); % Load #5
    AOP_file = uigetfile(file_spec, 'Choose AOP File'); % Load #8
    save(pressureFilenamesMat, 'LVP_file', 'AOP_file');
else
    load(pressureFilenamesMat);
end

% LVP 
data = load(sprintf('%s/%s', directory, LVP_file)); 
LV_ECG = data(:,4);
LV_Pressure = data(:,5);
LV_samplerate = 240;
LV_time = [1:1:length(LV_Pressure)]'/SAMPLE_RATE;

% AOP 
data = load(sprintf('%s/%s', directory, AOP_file)); 
AO_ECG = data(:,4);
AO_Pressure = data(:,5);
AO_samplerate = 240;
AO_time = [1:1:length(AO_Pressure)]'/SAMPLE_RATE;


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
    
%     % Get cut off frequencies for filtering
%     prompt = {'Enter cut-off frequency for low-pass filter: LVP', 'Enter cut-off frequency for low-pass filter: AOP'};
%     dlg_title = 'Input';
%     num_lines = 1;
%     defaultans = {'7.5', '6.25'};
%     cut_off_frequency = inputdlg(prompt, dlg_title, num_lines, defaultans);
    
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
    title('Filtered Aortic Pressure');
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
[LV_R_positions, LV_R_times, ~] = find_RR(LV_time, LV_ECG, LV_samplerate );
% AOP
[AO_R_positions, AO_R_times, ~] = find_RR(AO_time, AO_ECG, AO_samplerate );

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
        
        % If a minimum (eIVC) was found
        if ~isempty(l)
            eIVC(i) = l(1); % Get first local minimum
        else % Else, save as 1 and get rid of that AOP cycle afterwards
            eIVC(i) = 1;
        end
        
        % ES = dicrotic notch in aortic pressure trace
        %%%%% REQUIRES USER INPUT %%%%%
        ES(i) = getES_dicroticNotch(AO_cycles{2,i});
        
        % Get times corresponding to eIVC and ES
        eIVC_t(i) = AO_cycles{1,i}(eIVC(i));
        ES_t(i) = AO_cycles{1,i}(ES(i));
        
    end
    
    % Remove AOP cycles in which eIVC could not be found
    AO_cycles(:,(eIVC == 1)) = [];
    
    
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
%[RR_mean, TS_mean, TS_std, numFrames] = ExtractTriggerTime(model_directory);

% Get RR interval from CIM model file - IF NO TRIGGER TIMES ARE AVAILABLE IN HEADER, USE THIS LINE INSTEAD :-)
[RR_mean, numFrames] = getRR_CIM(model_directory);
if RR_mean > 0
    mri_time = linspace(1,RR_mean,numFrames)/1000;
else
    mri_time = linspace(1,1000,numFrames)/1000;
end

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

timepoints_MRI = [edMRI, eivcMRI, esMRI, eivrMRI, dsMRI, numFrames];

% If all five timepoints are available
if ~isnan(timepoints_MRI)
%if ~any(timepoints_MRI < 1)
    
    [LVP_average_unshifted, MRI_LVP_ed2eivc, MRI_LVP_eivc2es, MRI_LVP_es2eivr, MRI_LVP_eivr2ds, MRI_LVP_ds2ed] = ...
        InterpolatePressure_5Points(timepoints_MRI, LV_cycles, AO_cycles, ED, eIVC, ES, eIVR, DS, maxLVP);

% If only three timepoints are available: ED, ES, and DS
else
%elseif length(timepoints_MRI(timepoints_MRI==0)) == 2
    
    % Just keep points which are not NaN's
    %timepoints_MRI(timepoints_MRI==NaN) = [];
    timepoints_MRI(isnan(timepoints_MRI)) = [];
    
    [LVP_average_unshifted, MRI_LVP_ed2es, MRI_LVP_es2ds, MRI_LVP_ds2ed] = InterpolatePressure_3Points(timepoints_MRI, LV_cycles, AO_cycles, ...
        ED, ES, DS, maxLVP);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get final pressure trace to use as boundary conditions for parameter estimation

% Shift pressure so that LVP(DS) = 0
DS_pressure = min(LVP_average_unshifted); % Diastasis pressure
LVP_average = LVP_average_unshifted - DS_pressure; % Shifted LV pressure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create plots

% If all five timepoints are available
if length(timepoints_MRI) == 6
    
    plots_5points(timepoints_MRI, mri_time, LVP_average_unshifted, MRI_LVP_ed2eivc, ...
        MRI_LVP_eivc2es, MRI_LVP_es2eivr, MRI_LVP_eivr2ds, MRI_LVP_ds2ed, V, directory);

% If only three timepoints are available: ED, ES, and DS
elseif length(timepoints_MRI) == 4
    
    plots_3points(timepoints_MRI, mri_time, LVP_average_unshifted, MRI_LVP_ed2es, ...
        MRI_LVP_es2ds, MRI_LVP_ds2ed, V, directory);
    
end



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
