
close all

N = 64; % Frame size
L = 100; % Resampling
V_SOUND = 340.0;

% Mic positions
d = 0.024; % mic spacing in equilateral triangle
r_12 = [d 0 0]';
r_13 = [d/2 d*sqrt(3)/2 0]';

% Load signal
[y(:,1), Fs] = audioread('mic1.wav');
[y(:,2), Fs] = audioread('mic2.wav');
[y(:,3), Fs] = audioread('mic3.wav');
n = 0:length(y)-1;

% Storage for plots
main_storage = zeros(size(y));
high_storage = zeros(size(y));
detection = zeros(length(y),3);

% Filters
[b_main, a_main] = butter(4, 10000/(0.5*Fs),'high');
zi_main = [];

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

    % Calculate levels
    tap = 0;
    main_dB = db(rms(frame_main));
    if main_dB - mainAvg_dB > minMain_dB
        detection(i:i+N-1,2) = 0.8;
        tap = 1;
    end
    if main_dB - prevMain_dB > minMainIncrease_dB
        detection(i:i+N-1,1) = 0.5;
        tap = tap + 1;
    end
    
    prevMain_dB = main_dB;
    
    if tap == 2
        detection(i:i+N-1,3) = 0.4;
        
        % Tap detected
        
        % Resample and cross correlate
        frame_resampled = interpft(frame_main,N*L);
        maxlag = floor((0.024/340)*Fs*L);
        [r, lags] = xcorr(frame_resampled, maxlag);
        [maxcross, at_index] = max(r);
        
        % Extract delays (in meters...)
        d(1,1) = -lags(at_index(2))*V_SOUND/(Fs*L); % delay from 1 to 2
        d(2,1) = lags(at_index(7))*V_SOUND/(Fs*L); % delay from 1 to 3
        d(3,1) = -lags(at_index(6))*V_SOUND/(Fs*L); % delay from 2 to 3
        
        % Calculate normal vector of wavefront approximated by a plane
        C = [    r_12'    ;
                 r_13'    ;
             (r_13-r_12)' ];
        
        A = [0 0 1];
        b = -1;
        [v,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b);
        v = v/norm(v)
        
        
        figure
        subplot(3,1,1)
        plot(lags,r(:,2));
        title('12');
        subplot(3,1,2)
        plot(lags,r(:,6));
        title('23');
        subplot(3,1,3)
        plot(lags,r(:,7));
        title('31');
        pause
    end

    
    
    
    
    
    
    
    % Update moving averages
    mainAvg_dB = main_dB*mainAvg_alpha + mainAvg_dB*(1-mainAvg_alpha);
    
    
    % Store for plotting
    main_storage(i:i+N-1,:) = frame_main;

    
end


figure
h = zoom();
h.Motion = 'horizontal';
hold on


% Original
ax(1) = subplot(4,1,1);
plot(n,y);
hold on
plot(n,detection(:,3),'r');
title('Original');


% Main tap
ax(2) = subplot(4,1,2);
plot(n,main_storage);
title('Main');

% Detection
ax(3) = subplot(4,1,3);
plot(n,detection(:,1));
title('Increase');

ax(4) = subplot(4,1,4);
plot(n,detection(:,2));
title('Absolute');






% Link time axes
linkaxes(ax, 'x');      % Link all axes in x










