close all
tic

Fs = 48000;
N = 32; % Frame size
L = 1000; % Resampling
V_SOUND = 340.29;


H = dsp.AudioRecorder;
H.DeviceName = 'US-1800';
H.SampleRate = Fs;
H.SamplesPerFrame = N;
H.NumChannels = 6;
H.DeviceDataType = '24-bit integer';
H.OutputNumOverrunSamples = true;

% Internal cluster design
d = 0.024; % mic spacing in equilateral triangle
r_12 = [d 0 0]';
r_13 = [d/2 d*sqrt(3)/2 0]';

% Cluster positioning
p_1 = []'; % Cluster 1 position
p_2 = []'; % Cluster 2 position

% Filters
[b_main, a_main] = butter(4, 10000/(0.5*Fs),'high');
zi_main = [];

% Thresholds
minMainIncrease_dB = 10;
minMain_dB = 20;

% Levels
prevMain_dB = [0 0 0 0 0 0];
mainAvg_dB = [0 0 0 0 0 0];
mainAvg_alpha = 0.010*N/64;

n = 0;
while true
    
    % Extract frame
    [frame , overrun] = step(H);
    if overrun > 0
       fprintf('Overrun\n'); 
    end

    % Main tap highpass filter
    [frame_main, zi_main] = filter(b_main,a_main,frame,zi_main);

    % Split into two clusters
    
    
    % Check levels for both clusters
    main_dB = db(rms(frame_main));
    for i = 1:3:4
        if main_dB(i:i+2) - mainAvg_dB(i:i+2) > minMain_dB
            if main_dB(i:i+2) - prevMain_dB(i:i+2) > minMainIncrease_dB

                % Tap detected
                toc
                fprintf('%i - Cluster %i: %f dB rms level, %f dB increase\n', n, i, mean(main_dB(i:i+2) - mainAvg_dB(i:i+2)), mean(main_dB(i:i+2) - prevMain_dB(i:i+2)));
                n = n+1;
                cluster_frame_main = frame_main(:, i:i+2);
                
                
                % Resample and cross correlate
                cluster_frame_resampled = interpft(cluster_frame_main,N*L);
                maxlag = floor((0.024/340)*Fs*L);
                [r, lags] = xcorr(cluster_frame_resampled, maxlag);
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
                [v,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,[],[],A,b);
                v = v/norm(v);

                % Figure out the x y coordinates, at 1m, for demo use
                screenPoint = v(1:2)/v(3);
                
                figure(i)
                scatter(screenPoint(1),screenPoint(2),'r');
                xlim([-1,1]);
                ylim([-1,1]);
                hold on
                drawnow
                
                % increase the avg, to avoid multiple detections
                mainAvg_dB(i:i+2) = main_dB(i:i+2);
                
                
            end
        end
    end
    prevMain_dB = main_dB;
    

    % Update moving averages
    mainAvg_dB = main_dB*mainAvg_alpha + mainAvg_dB*(1-mainAvg_alpha);
    
    
end

  
release(H);
    







