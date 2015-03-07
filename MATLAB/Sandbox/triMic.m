
close all

N = 64; % Frame size

% Load signal
[y(:,1), Fs] = audioread('mic1.wav');
[y(:,2), Fs] = audioread('mic2.wav');
[y(:,3), Fs] = audioread('mic3.wav');
n = 0:length(y)-1;

% Storage for plots
main_storage = zeros(size(y));
high_storage = zeros(size(y));
low_storage = zeros(size(y));
detection = zeros(length(y),2);

% Filters
[b_main, a_main] = butter(4, 10000/(0.5*Fs),'high');
zi_main = [];
[b_low, a_low] = butter(4, 8000/(0.5*Fs),'low');
zi_low = [];

% Thresholds
minMainIncrease_dB = 10;
minMain_dB = 20;

% Levels
prevMain_dB = 0;
mainAvg_dB = 0;
mainAvg_alpha = 0.015*N/64;

% Iterate over all frames
n_frames = floor(length(y)/N)-1;
for i = 1:N:n_frames*N
    
    % Extract frame from raw data
    frame = y(i:i+N-1,:);
    
    % Main tap highpass filter
    [frame_main, zi_main] = filter(b_main,a_main,frame,zi_main);

    % Pre tap lowpass filter
    %[frame_low, zi_low] = filter(b_low,a_low,frame,zi_low);

    
    % Calculate levels
    main_dB = db(rms(frame_main));
    if main_dB - prevMain_dB > minMainIncrease_dB
        detection(i:i+N-1,1) = 0.5;
    end
    if main_dB - mainAvg_dB > minMain_dB
        detection(i:i+N-1,2) = 0.8;
    end
    prevMain_dB = main_dB;
    
    % Update moving averages
    mainAvg_dB = main_dB*mainAvg_alpha + mainAvg_dB*(1-mainAvg_alpha);
    
    
    % Store for plotting
    main_storage(i:i+N-1,:) = frame_main;
    %low_storage(i:i+N-1,:) = frame_low;
    
    
end


figure
h = zoom();
h.Motion = 'horizontal';
hold on


% Original
ax(1) = subplot(4,1,1);
plot(n,y);
title('Original');

% Main tap
ax(2) = subplot(4,1,2);
plot(n,main_storage);
title('Main');

% Detection
ax(3) = subplot(4,1,3);
plot(n,detection(:,1));
title('Increase');

% Detection
ax(4) = subplot(4,1,4);
plot(n,detection(:,2));
title('Absolute');





% Link time axes
linkaxes(ax, 'x');      % Link all axes in x










