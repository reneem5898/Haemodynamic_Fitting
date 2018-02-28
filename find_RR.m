function [R_positions, R_times, RR] = find_RR(time, ECG, samplerate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses a modified algorithm of the version provided by
% librow.com to locate R peaks from an ECG trace
% 
% Adapted by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 11 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove lower frequencies
f_LV = fft(ECG);
f_LV(1:round(length(f_LV)*5/samplerate)) = 0;           % Removes lower end values, numerator originally*5
f_LV(end - round(length(f_LV)*5/samplerate):end) = 0;   % Removes higher end values
corrected = real(ifft(f_LV));

% Filter - first pass
WinSize = floor(samplerate * 571 / 1000); %% RM - why 571?
if rem(WinSize,2) == 0
    WinSize = WinSize+1;
end
filtered1 = ecgdemowinmax(corrected, WinSize);

% Scale filtered1
peaks1 = filtered1/(max(filtered1)/7); %% RM - why / 7???

% Filter by threshold filter
for i = 1:length(peaks1)
    if peaks1(i) < 4
        peaks1(i) = 0;
    else
        peaks1(i) = 1;
    end
end

% %ET40 - when negative t wave is the start of the ECG (9th data point), it is confused with R
% % peak - if a peak is detected within first 20 samples then ignore
% for i=1:20
%     if (peaks1(i)~=0)
%         peaks1(i)=0;
%     end
% end

R_positions = find(peaks1);
R_times = time(R_positions);

% %Check LV alignment
% figure
% subplot(2,1,1)
% plot(x_LV, corrected);
% title('LV ECG')
% subplot(2,1,2)
% plot(x_LV, peaks1);
% title('LV R peak times');

%Calculate RR intervals
for i = 1:(length(R_times)-1)
    RR(i) = R_times(i+1) - R_times(i);
end

end

