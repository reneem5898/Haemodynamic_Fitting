% Plot pressure and ECG from BioBeat data

% Directory
directory = 'N:\AMRG\AMRG\MRIresearch\kat\Pressure';

% Each dataset has 5 columns:
% Column 2 - ECG (mV)
% Column 5 - Pressure (mmHg)

%
% ADHB_ID	Number_Snaps	Study_ID	Pressure_Snap	BinaryFileName	OtherInfo	Comment
% CATH-3681	15		1	X000.txt	X000.inf	AO
% 			2	X001.txt	X001.inf	AO
% 			3	X002.txt	X002.inf	AO
% 			4	X003.txt	X003.inf	LV
% 			5	X004.txt	X004.inf	LV
% 			6	X005.txt	X005.inf	LV
% 			7	X006.txt	X006.inf	LV
% 			8	X007.txt	X007.inf	Pullback LV -> AO
% 			9	X008.txt	X008.inf	AO
% 			10	X009.txt	X009.inf	AO
% 			11	X010.txt	X010.inf	AO
% 			12	X011.txt	X011.inf	Arch
% 			13	X012.txt	X012.inf	ART
% 			14	X013.txt	X013.inf	ART
% 			15	X014.txt	X014.inf	ART
%

% Loop through and plot ECG and pressure for each dataset

for i = 1:14

    data = load(sprintf('%s/X%03d.txt', directory, i-1));
    
    ecg = data(:,2);
    pressure = data(:,5);
    
    FH = figure('position', [100 100 800 600]);
    subplot(2,1,1)
    plot(ecg)
    hold on
    title(sprintf('X%03d ECG',i-1), 'FontSize', 14)
    xlabel('Data Point', 'FontSize', 12)
    ylabel('Voltage (mV)', 'FontSize', 12)
    
    subplot(2,1,2)
    plot(pressure)
    hold on
    title(sprintf('X%03d Pressure',i-1), 'FontSize', 14)
    xlabel('Data Point', 'FontSize', 12)
    ylabel('Pressure (mmHg)', 'FontSize', 12)

    saveas(FH, sprintf('%s/X%03d.png', directory, i-1));
    close(FH);
end

% Test peak finding from ECG
samplerate = 240;