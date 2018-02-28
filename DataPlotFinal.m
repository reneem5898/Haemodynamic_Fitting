%This file 
% 1) determines the positions of the R peaks of LV ECG data from
% OHIO University
% 2) low pass filters pressure data from OHIO University after analysing
% the power spectrum
% 3) extracts all individual cardiac cycle using the R-R interval
% 4) calculates a representative pressure trace by averaging nearest pressure data
% from all individual cycles using a time scale defined by the average RR
% interval and spaced at 1 ms (sample rate of data)
% 5) calculates first and second derivative of LVP using forward or central
% difference
% 6) passes the second derivative through a 5 point median filter to
% remove 'salt and pepper' noise
% 7) identifies cardiac events: endIVR, ED and DS using LVP data
% 8) extracts pressure data using piecewise linear scaling between marked
% cardiac events to plot PV loop

clear
clf
close all

%% Extract data
%%----USER PARAMETER: file to read-----
%file=input('Please enter name: ', 's');
%filetype=input('Please enter filetype (1 for ET40,41, 2 for ET28,ET46): ', 's'); 
file = 'ET40MI';                      
filetype = '1';
deriv_method = 'central';                                               
                                    
[ period_LV, x_LV, LV_Pressure, LV_ECG] = load_data( file, filetype );
samplerate = 1/period_LV;

%% Find R peaks
[LV_R_positions, LV_R_times, RR_LV] = find_RR( x_LV, LV_ECG, samplerate );

plot(x_LV, LV_Pressure, '*', 'markersize', 2);
title('Raw LV Pressure');
ylabel('Pressure (mmHg)');
xlabel('Time (ms)');

%%Power Spectrum
PS_LVP = LV_Pressure; 
PS_xLV = x_LV;

%Fourier transform
L = length(PS_LVP);
NFFT = 2^nextpow2(L);
PS = fft(PS_LVP, NFFT)/L;
x_PS = samplerate(1)/2*linspace(0,1,NFFT/2+1);

%Plot power spectrum
figure
plot(x_PS, 2*abs(PS(1:NFFT/2+1)));
title('Power Spectrum')
ylabel('Amplitude')
xlabel('Frequency (Hz)')

%% Filtering
[LV_Pressure] = filter_pressure(LV_Pressure, samplerate);

%% Extract all cardiac cycles
[ LV_cycles] = extract_individual_cycles( x_LV, LV_R_positions, RR_LV, LV_Pressure, LV_ECG);

% Check: Overlay all LV pressure traces
figure
for i = 1:length(RR_LV)
    plot(LV_cycles{1,i},LV_cycles{2,i}, '*', 'markersize', 2);
    hold on;
end
title('All LV pressure traces')
hold off;

%% Representative Trace
[ representative_LV_times, representative_LV_pressure, representative_LV_ECG] = find_representative_trace(RR_LV, LV_cycles, samplerate);

%% Smooth unscaled LVP to mark cardiac events
h=representative_LV_times(2)-representative_LV_times(1);              %step size for calculating derivatives

[deriv1_LV, deriv2_LV]=find_derivative(h, representative_LV_pressure, deriv_method); %find derivatives
smooth_deriv2LV=medfilt1(deriv2_LV, 5); %3 point not good enough for ET40

% Mark cardiac event                                                    
[diastasis_index, endIVR_index, ED_index ] = mark_cardiac_events(representative_LV_pressure, deriv1_LV, smooth_deriv2LV); 

%Plot derivatives
% figure 
% if strcmp(deriv_method, 'forward')
%     subplot(4,1,1);
%     plot(representative_LV_times, representative_LV_pressure, '*');
%     subplot(4,1,2);
%     plot(representative_LV_times(1:end-1), deriv1_LV, '*');
%     subplot(4,1,3);
%     plot(representative_LV_times(1:end-2), deriv2_LV, '*');
%     subplot(4,1,4);
%     plot(representative_LV_times(1:end-2), smooth_deriv2LV, '*');
% elseif strcmp(deriv_method, 'central')
%     subplot(4,1,1);
%     plot(representative_LV_times, representative_LV_pressure, '*');
%     subplot(4,1,2);
%     plot(representative_LV_times, deriv1_LV, '*');
%     subplot(4,1,3);
%     plot(representative_LV_times, deriv2_LV, '*'); 
%     subplot(4,1,4);
%     plot(representative_LV_times, smooth_deriv2LV, '*'); 
% end

diastasis_time = representative_LV_times(diastasis_index);
endIVR_time = representative_LV_times(endIVR_index);
ED_time = representative_LV_times(ED_index);

%% Plot results
        
% Representative pressure and ECG with marked cardiac event
figure
subplot(2,1,1)
plot(representative_LV_times, representative_LV_pressure, '*', 'markersize', 2);
hold on;
line([endIVR_time endIVR_time], [-50 100],'Color',[1 0 0])
line([diastasis_time diastasis_time], [-50 100], 'Color', [0 1 0])
line([ED_time ED_time], [-50 100], 'Color', [0 0 0]);
title('LV Pressure Trace with Cardiac Events');
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
legend('LV Pressure Trace','End IVR','Diastasis','End Diastole', 'Location', 'southwest')
hold off;
subplot(2,1,2);
plot(representative_LV_times, representative_LV_ECG, '*', 'markersize', 2);
title('LV Representative ECG');
ylabel('Voltage (mV)')
xlabel('Time (s)')

% Representative pressure, 1st derivative, 2nd derivative and ECG with
% marked events
figure 
if strcmp(deriv_method, 'forward')
    subplot(4,1,1);
    plot(representative_LV_times, representative_LV_pressure, '*', 'markersize', 2);
    title('Representative LVP with marked cardiac events');
    ylabel('LVP (mmHg)')
    xlabel('Time (s)')
    subplot(4,1,2);
    plot(representative_LV_times(1:end-1), deriv1_LV, '*', 'markersize', 2);
    title('First Derivative');
    ylabel('mmHg/s')
    xlabel('Time (ms)')
    subplot(4,1,3);
    plot(representative_LV_times(1:end-2), smooth_deriv2LV, '*', 'markersize', 2);
    title('Smoothed Second Derivative');
    ylabel('mmHg/s^2')
    xlabel('Time (s)')
    subplot(4,1,4);
    plot(representative_LV_times, representative_LV_ECG);
    title('LV Representative ECG');
    ylabel('Voltage (mV)')
    xlabel('Time (s)')
elseif strcmp(deriv_method, 'central')
    subplot(4,1,1);
    plot(representative_LV_times, representative_LV_pressure, '*', 'markersize', 2);
    title('Representative LVP with marked cardiac events');
    ylabel('LVP (mmHg)')
    xlabel('Time (s)')
    subplot(4,1,2);
    plot(representative_LV_times, deriv1_LV, '*', 'markersize', 2);
    title('First Derivative');
    ylabel('mmHg/s')
    xlabel('Time (ms)')
    subplot(4,1,3);
    plot(representative_LV_times, smooth_deriv2LV, '*', 'markersize', 2); 
    title('Smoothed Second Derivative');
    ylabel('mmHg/s^2')
    xlabel('Time (s)')
    subplot(4,1,4);
    plot(representative_LV_times, representative_LV_ECG);
    title('LV Representative ECG');
    ylabel('Voltage (mV)')
    xlabel('Time (s)')
end
for i=1:4
   subplot(4,1,i);
   hold on;
   yL = get(gca,'YLim');
   line([endIVR_time endIVR_time], yL,'Color',[1 0 0])
   line([diastasis_time diastasis_time], yL, 'Color', [0 1 0])
   line([ED_time ED_time], yL, 'Color', [0 0 0]);
   legend('data', 'endIVR', 'DS', 'ED', 'Location', 'BestOutside');
   hold off;
end

%% Extract Pressure information for PV loop
events = 0;
totalFrames=input('Please enter total MRI frames: ');
eventIndices=[];            %stores time index of cardiac events for pressure trace
MRIFrames=[];               %stores frames of cardiac events used
while events < 2
    DS=input('Do you want to use Diastasis to fix PV loop? 1 for yes, 0 for no: ');
    if (DS)
        eventIndices=[diastasis_index];
        events=events+1;
        DS_MRI=input('Please specify MRI frame for when diastasis was observed: ');
        MRIFrames=[DS_MRI];
    end
    ED=input('Do you want to use end diastole to fix PV loop? 1 for yes, 0 for no: ');
    if (ED)
        eventIndices = [eventIndices ED_index];
        events = events+1;
        ED_MRI = input('Please specify MRI frame for when ED was observed: ');
        MRIFrames = [MRIFrames ED_MRI];
    end
    endIVR = input('Do you want to use endIVR to fix PV loop? 1 for yes, 0 for no: ');
    if (endIVR)
        eventIndices = [eventIndices endIVR_index];
        endIVR_MRI = input('Please specify MRI frame for when endIVR was observed: ');
        MRIFrames = [MRIFrames endIVR_MRI];
        events = events+1;
    end
    initial=input('Do you want to use initial MRI frame to fix PV loop? 1 for yes, 0 for no: ');
    if (initial)
        initialframeIndex = input('Please specify total delay time after R wave for when first MRI frame was taken: ');
        eventIndices = [eventIndices initialframeIndex];
        MRIFrames = [MRIFrames 1];
        events = events+1;
    end
end

% Combine pressure indices with MRI frames
PV_eventIndices(1,:) = eventIndices;
PV_eventIndices(2,:) = MRIFrames;

% Sort pressure indices in ascending order with correspondence
[~,p] = sort(PV_eventIndices(1,:),2);
PV_eventIndices = PV_eventIndices(:,p);

PV_indices = [];          %%pressure index for PV loop for each MRI frame
% Linearly align pressure and volume between events used (need two or more
% events)
for i = 1:(events)
    if (i == events)
        framesBetween = PV_eventIndices(2, 1)+totalFrames-PV_eventIndices(2,end);
        temp_index = PV_eventIndices(1,1)+length(representative_LV_pressure);
        last_pressure_indices = (round(linspace(PV_eventIndices(1, end), temp_index, framesBetween+1)));
        for j = 1:length(last_pressure_indices)
           while (last_pressure_indices(j) > length(representative_LV_pressure)) 
               last_pressure_indices(j) = last_pressure_indices(j)-length(representative_LV_pressure);
           end
        end
        PV_indices = [PV_indices(1:end) last_pressure_indices(2:end-1)];
    else
        framesBetween = PV_eventIndices(2, i+1)-PV_eventIndices(2, i);
        PV_indices = [PV_indices(1:end-1) round(linspace(PV_eventIndices(1, i), PV_eventIndices(1, i+1), framesBetween+1))];
    end
end

PV_pressure = representative_LV_pressure(PV_indices);
PV_pressure = PV_pressure - representative_LV_pressure(diastasis_index);    % offset so diastasis is at pressure=0

% Representative LV Pressure with markers where pressure data has been extracted
figure
plot(representative_LV_times, representative_LV_pressure, '*', 'markersize', 2);
hold on;
line([endIVR_time endIVR_time], [-50 100],'Color',[1 0 0])
line([diastasis_time diastasis_time], [-50 100], 'Color', [0 1 0])
line([ED_time ED_time], [-50 100], 'Color', [0 0 0]);
plot((PV_indices-1)/1000, representative_LV_pressure(PV_indices), '*', 'markersize', 5, 'Color', [1 0 0]);
axis([0 0.7 -20 120]);
title('LV Pressure Trace with Cardiac Events');
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
legend('LV Pressure Trace','End IVR','Diastasis','End Diastole', 'Location', 'BestOutside')
