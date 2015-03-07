function [y, zf] = tapFilter(x,F0,Fs,BW,N,zi)
%FILTER Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 8.4 and the DSP System Toolbox 8.7.
% Generated on: 11-Feb-2015 19:09:53



h = fdesign.peak('N,F0,BW', N, F0, BW, Fs);

Hd = design(h, 'butter', ...
    'SOSScaleNorm', 'Linf');



[y, zf] = filter2(Hd,x,zi);

