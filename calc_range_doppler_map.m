function [numRangeBins,numDopplerBins,doppler_FFT] = calc_range_doppler_map(c,N,Zeropad,T,BW,kf,c0,oRS)
c = c(:,1,:);
%% Maximum Range calculation settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% r = (c0*T*oRS.fSamplingRate)/2*BW
% range-Doppler map
repeat=double(oRS.oEPRadarBase.num_chirps_per_frame);                                     
yscale=2*kf/c0;                                                                           % scaling for the y-axis as a distance in m   
ydata=linspace(0,1-1/Zeropad,Zeropad)*double(oRS.oEPRadarSoli.samplerate_Hz)/yscale;      % 0-3.2m
yaxis=linspace(0,ydata(end)/2,Zeropad/2);                                                 % scaling from 0-1.6 with 256 evenly spaced points   



%% Maximum Velocity calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vmax = (lambda)/(4*chirpDuration)
% range-doppler map
fc=(double(oRS.oEPRadarFMCW.lower_frequency_kHz+oRS.oEPRadarFMCW.upper_frequency_kHz)/2)*1e3; % center frequency
lambda=c0/fc;                                                                                 % = 0.005 radar wavelength in m
max_velocity=lambda/(4*T);                                                                    % = 19.5313
xaxis=linspace(-max_velocity,max_velocity,Zeropad/2);                                         % scaling from -19.5131 to 19.5131 with 256 evenly spaced points 


%% Perform range and Doppler FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% common
% hTime=figure;
Hann_window=hann(N,'periodic');                                                 % hann(Length of window, 'periodic');
ScaleHannWin = 1/sum(Hann_window);                                              % replaces the scaling of the fft with length of spectrum with sum(window) 

Cheb_window=chebwin(repeat,150);                                                % chebwin(Length of window, difference between sidelobe mag to mainlobe mag)
ScaleChebWin=1/sum(Cheb_window);

%%  º∆À„range-Doppler map
for count=1:repeat
range_FFT(count,:)=(fft(c(:,:,count)'.*Hann_window',Zeropad).*ScaleHannWin);    % For-loop for calculation of chirp 1 to 40
end
%Coherent integration and subtraction 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mean_row=mean(range_FFT);                                                  
for kk=1:Zeropad                                                                % from 1 to 512
    range_FFT(:,kk)=range_FFT(:,kk)-Mean_row(kk);                               % DC-Removal  offset at 0Hz:                                                                         
end
for jj=1:Zeropad/2                                                              % from 1 to 256                                                            
     doppler_FFT(jj,:)=db(abs(fftshift(fft(range_FFT(:,jj)'.*Cheb_window',Zeropad/2).*ScaleChebWin)));  
end
numRangeBins = length(yaxis);
numDopplerBins = length(xaxis);
end