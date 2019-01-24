function [rvRange, rvAngle, rvMax] = AoADetection(mxData, scF0, scF1, scTsw, scSpace, scNfft, rvBP, scC0)
% function [rvRange, rvAngle, rvMax] = AoADetection(mxData, scF0, scF1, scTsw, scSpace, scNfft, scC0)
% mxData: raw data input from a uniform circular array or uniform linear array, one sweep per column
% scF0: start frequency in Hz (e.g. 57e9)
% scF1: stop frequency in Hz (e.g. 57e9+5.8e6 = 57,005,800,000)
% scTsw: sweep time in s (e.g. 150e-6)
% rvBP: normalized cutoff frequency for bandpass filter from 0...1 with 1^= nyquist rate梙alf the sample rate
% scSpace: spacing between two receivers in m (e.g. scWL*.7 for scWL=5e-3, F=6GHz)
% scNfft: fft length with zero padding (e.g. 2^12 for 256 sweep taps)
% scC0: speed of light in m/s (i.e. 3e8)
%
% rvRange: % row vector with detected range values in m  
% rvAngle: % row vector with detected angle of arrival (AoA) in rad
% rvMax: % normalized maximum value of peak, 1 ^= full scale of analysed sine
%

% compute physical constants and hardware properties
scN = size(mxData,1); % no. of sweep taps (256)

if nargin<8
    scC0 = 3e8; % speed of light in m/s
end
if nargin<7
    rvBP = [0 1]; % no band pass
end
if nargin<6
    scNfft = scN; % no zero padding
end

% hardware settings
scB = scF1-scF0; % sweep bandwidth in Hz (5.8e6)
scWL = scC0/scF1; % maximum signal wavelength in m

scTs = scTsw/scN; % sample duration in s, i.e. scN samples in one sweep interval
scSR = 1/scTs; % sample rate in 1/s

% FFT
mxData = mxData.*repmat(hann(size(mxData,1)),1,4);
mxFFT = fft(mxData,scNfft,1); % fft with zero padding across dimension 1, i.e. taps/time

mxFFTHalf =(mxFFT(1:scNfft/2,:).*repmat((0:scNfft/2-1)',1,size(mxFFT,2))); % take only half spectrum and compensate free field damping

% apply band pass
loCutInd = round(rvBP(1)*scNfft/2)+1; % index based, i.e. +1 50
hiCutInd = round(rvBP(2)*scNfft/2)+1; % index based, i.e. +1 417
mxFFTHalf(1:loCutInd,:) = mxFFTHalf(1:loCutInd,:).*repmat(linspace(0,1,loCutInd)',1,4);
mxFFTHalf(hiCutInd:end,:) = 0;

% range and angle of arrival detection
% 得到最大峰点值的index和相位差
[rvPeaks, rvPhDelt] = FFT2PeakNPhase(mxFFTHalf); % find peaks and phase differences


rvRange = getRange(rvPeaks, scNfft, scSR, scB/scTsw, scC0); % detected range value in m
% disp(rvRange);
rvAngle = getAoA(rvPhDelt, scSpace/scWL); % detected AoA in rad
rvMax = 2.*mxFFTHalf(rvPeaks+1)./scN;
% 去除noise
if max(abs(rvMax))<2 % threshold for noise floor
    rvRange = zeros(1,4);
    rvAngle = zeros(1,4);
end    
end


function [ rvPeaks, rvPhDelt ] = FFT2PeakNPhase( mxFFT )
% mxFFT: half spectrum of an fft
%
% rvPeaks: zero based bin index of peak position
% rvPhDelt: phase information of peak position
%

% find peaks in spectrum
% [rvPks,rvLocs] = findpeaks(abs(mxFFT(:,1))'); % get all peaks (req. Signal Processing Toolbox)
% 得到峰值点index
[~,rvPeakInd] = max(abs(mxFFT)); % get max peak
rvPeaks = rvPeakInd-1; % zero based

% phase difference
rvLinInd = zeros(2,size(mxFFT,2));
% 取出峰点值
rvLinInd(1,:) = sub2ind(size(mxFFT), rvPeakInd, 1:size(mxFFT,2)); % linear index to pick out peaks from matrix
rvLinInd(2,:) = sub2ind(size(mxFFT), rvPeakInd, [2:size(mxFFT,2),1]); % linear index to pick out peaks from matrix
% 去除重复周期，如2pi
rvPh = unwrap(angle(mxFFT(rvLinInd)),[],1); % get phase information only
% 得到相位差
rvPhDelt = diff(rvPh,1,1); % get phase deltas
end


function cvRange = getRange(cvBins, scNfft, scSR, scK, scC0)
% cvBins: FFT bins starting with 0
% scNfft: FFT size in Taps
% scSR: sample rate in 1/s
% scK: 频率带宽/扫描周期 = 频率变化率
% scK: steepness, i.e. = B/Tsw in Hz/s with B: bandwidth and Tsw: sweep duration
% scC0: signal speed in m/s, default: scC0 = 3e8; % speed of light
%
% cvRange: distance in m
%

if nargin<5
    scC0 = 3e8; % speed of light in m/s
end

cvPeakF = (cvBins).*(scSR/scNfft); % peak frequency in Hz
cvPeakT = cvPeakF./scK; % peak round trip time in s
% peak一个来回的时间的1/2*波速=range value
cvRange = (cvPeakT./2).*scC0; % detected range value in m
end


function cvAlph = getAoA(cvPhDelt, scXi)
% cvPhDelt: phase deltas in rad
% scXi = scSpace(两个receivers之间的距离)/scWL(最大信号波长)
% scXi: normalized space between receivers, i.e. = S/WL with S: space in m and WL: wavelength in m 
%
% cvAlph: AoA in rad, when looking from the receiver to the source:
% 0: broadside, -pi/2: left, pi/2: right
%

cvAlphSin = (cvPhDelt/(2*pi*scXi)); % normalize to xi (scSpace/scWL)
% asin-- 反正弦函数，返回弧度
cvAlph = asin(cvAlphSin);
end