close all

Fs = 48000;
N = 64; % Frame size
L = 100; % Resampling
V_SOUND = 340.0;


H = dsp.AudioRecorder;
H.DeviceName = 'US-1800';
H.SampleRate = Fs;
H.SamplesPerFrame = 64;
H.NumChannels = 3;
H.DeviceDataType = '24-bit integer';
H.OutputNumOverrunSamples = true;



% Extract frame from raw data
[frame , overrun] = step(H);



release(H);