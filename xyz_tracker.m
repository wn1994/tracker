dbstop if error
clc
disp('******************************************************************');
addpath('..\..\RadarSystemImplementation');                           % add MATLAB API
clear all
close all                                                                      % close and delete ports
warning off
%% setup object and show properties
szPort = findRSPort;                                                           % change the available port
oRS = RadarSystem(szPort);                                                     % setup object and connect to board
disp('oRS object - properties before set block:');
oRS

%% set properties
% Changing some default properties on the board for proper operation of the code.
oRS.oEPRadarFMCW.lower_frequency_kHz = oRS.oEPRadarBase.min_rf_frequency_kHz+1e6; % = 57 GHz
oRS.oEPRadarFMCW.upper_frequency_kHz = oRS.oEPRadarBase.min_rf_frequency_kHz+7e6; % = 63 GHz
oRS.oEPRadarFMCW.direction = 0;                                                   % 'up-chirp'
oRS.oEPRadarFMCW.tx_power = 31;                                                   % valid powerrange: 0 - 31
% 10 will set the conducted transmit power to a value of 0 dBm
oRS.oEPRadarSoli.samplerate_Hz = 2000000;                                         % 2 MHz max. sampling rate
oRS.oEPRadarBase.num_chirps_per_frame = 31;
% oRS.oEPRadarBase.num_chirps_per_frame = 1;                                        % one chirp per radar transmission
oRS.oEPRadarBase.num_samples_per_chirp = 128;                                     % 128 sample points per chirp
oRS.oEPRadarBase.rx_mask = bin2dec('1100');                                       % turn on up to 4 receiverchannels
% every bit represents one receiver, bit mask is oRS.sRXMask=[RX4 RX3 RX2 RX1]
oRS.oEPRadarSoli.tx_mode = 0;                                                     % TX1 transmitter only
% oRS.oEPRadarSoli.tx_mode = 2;                                                   % 2 for two transmitter sequencially
disp('oRS object - properties after set block:');


N=double(oRS.oEPRadarBase.num_samples_per_chirp);
Zeropad=N*4;                                                                   % = Zeropadding to multiple of N where N*4 is good choice
T=double(oRS.oEPRadarBase.chirp_duration_ns)*1e-9;                             % chirp time in seconds
BW=double(oRS.oEPRadarFMCW.upper_frequency_kHz-oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % Radar bandwidth in Hz
kf=BW/T;                                                                       % slope rate in Hz/s;
c0=3e8;                                                                        % speed of light im m/s;

scWL=c0/(double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3);                    % maximum wavelength
scSpace=0.7*scWL;                                                              % spacing between the receivers
xscale=2*kf/c0;                                                                % scaling for the x-axis as a distance in m
fadc =double(oRS.oEPRadarSoli.samplerate_Hz);
xdata=linspace(0,1-1/Zeropad,Zeropad)*fadc/xscale;                             % scaling the x-axis to wproper range to frequency
xdata1=linspace(0,1-1/N,N)*fadc/xscale;                                        % scaling for

Hann_window=hann(N,'periodic');                                                % hann(Length of window, 'periodic');
ScaleHannWin = 1/sum(Hann_window);                                             % replaces the scaling of the fft with length of spectrum with sum(window)


%% start Angle of Arrival estimation
% starting the while-loop where the process of AoA is computed and collecting raw_data from the board
% �ϴα����rangeֵ
lastRange = 0;
lastAoA = 0;
x=0;
y=0;
xx=0;
yy=0;
aa = 0;
bb = 0;
% ȥ��ǰһС������
cutLength = 10;

dataSize = 0;
bufferSize = 2000;
% red ����
mxHistory = NaN(2,bufferSize);
% blue ����
mxAllHistory = NaN(2,bufferSize);
% �洢range��AOA
RangeAndAoAdata = NaN(2,bufferSize);
hTime = figure;
set(hTime,'units','normalized','position',[0.3 0.35 0.35 0.50]);
% ���귶Χ
limValue = 60;
while ishandle(hTime)
    pause(0.000001); %����Ҫ������� Ҫ��Ȼ��������޷��õ���������
    if isletter(get(hTime,'CurrentCharacter'))
        break;
    end
    % trigger radar chirp and get radar raw data
    [c, ~]=oRS.oEPRadarBase.get_frame_data;
    % c = [cc(:,1,1), cc(:,2,1)];
    
    % DC Removal
    % accomplished by substracting the mean of the data from the raw data suppressing the DC signal and filtering the data with a high-pass to remove all low frequencies
    % ͨ��������DC�źŵ�ԭʼ�����м�ȥ���ݵ�ƽ��ֵ���ø�ͨ�˲��������������е�Ƶ
    
    c(:,1)=c(:,1)-mean(c(:,1));
    c(:,2)=c(:,2)-mean(c(:,2));
    
    Nbut=4;                                                                        % 4th order of filter
    Wn=.05;                                                                        % cutoff frequency
    [b,a]=butter(Nbut,Wn,'high');                                                  % butter(order of filter, cutoff frequency, type of filter);
    % ��ͨ�˲�
    data1=filter(b,a,c(:,1)');                                                     % filtering
    data2=filter(b,a,c(:,2)');
    
    
    %% Prepare AoA
    % Computing the fast-fourier-transform of the raw data and finding any peaks for a region of data,
    %here the data has to be in a range of 0.1 to 0.9. Peaks are than used for determing range and signal strength of object
    
    FFT_dat1=fft(data1.*Hann_window',Zeropad);
    scaled_FFT1=db(abs(FFT_dat1(1:length(xdata))*ScaleHannWin));                   % FFT
    
    FFT_dat2=fft(data2.*Hann_window',Zeropad);
    scaled_FFT2=db(abs(FFT_dat2(1:length(xdata))*ScaleHannWin));
    IF_info=[FFT_dat1' FFT_dat2'];                                                % Create a matrix with two receivers
    xmin = 0.1;                                                                % Minimum range value
    xmax = 0.5;                                                                % Maximum range value
    region_of_interest = xmax>xdata & xdata>xmin;                              % Desired range value for peak detection, output logic array either '0' or '1' each cell
    start_bin=find(region_of_interest,1)-1;                                    % find first cell with '1' and save cellnumber - 1
    
    [rvPks,rvLocs] = findpeaks(scaled_FFT1(region_of_interest),'MINPEAKHEIGHT',-35,'MINPEAKDISTANCE',5);
    % find local peaks only for cells with logic '1'
    
    iter=length(rvLocs);
    
    
    if(iter>0)
        dataSize = dataSize + 1
        figure(hTime);
        clf
        hold on
        for kk=1:iter
            Desired_bin(kk) = start_bin+rvLocs(kk);                                 % Desired bin
            target_range(kk)= xdata(Desired_bin(kk));                               % Target range value
            target_signal_value(kk)=scaled_FFT1(Desired_bin(kk));                   % Target signal level
            
            
            %%  Angle of Arrival estimation
            % AoA is determined by taking the phased difference between the two receivers and calculating the AoA:
            %  ͨ��ȡ����������֮������ȷ��AoA��������AoA
            % AoA = arcsin((delta_phi/2pi)*(max. wavelength/spacing))
            
            rvPh(kk,:) = unwrap(angle(IF_info(Desired_bin(kk),:)),[]);                    % get phase information from the desired target bin for both receivers
            rvPhDelt(kk) = diff(rvPh(kk,:));                                              % get phase difference between the two receivers
            scAlphSin(kk) = (rvPhDelt(kk)/(2*pi))*(scWL/scSpace);                         % Angle of arrival formula alpha = (delta_phase/2*pi)*(lambda/0.7*lambda)
            Angle_of_arrival(kk) = asind(scAlphSin(kk));                                  % Angle value in degrees
        end
        
        
        %% plot tracker
        % �ų��쳣��
        % �ҳ���⵽�ĵ�������һ�����������ĵ��index
        [~,lastRangeIdx] = min(abs(target_range - lastRange));
        
        %% �������м�⵽�Ĺ켣����ɫ��
        if dataSize >= cutLength
            xx=100 * target_range(1) * sind(Angle_of_arrival(1));
            yy=100 * target_range(1) * cosd(Angle_of_arrival(1));
            mxAllHistory(:,dataSize) = [xx;yy];
            plot(xx, yy, 'or', 'MarkerFaceColor', 'blue', 'MarkerSize', 15);
            hold on;
            plot(mxAllHistory(1,1:dataSize), mxAllHistory(2,1:dataSize), 'blue', 'LineWidth', 2);
            hold on;
        end
        
        %% �������˺�Ĺ켣����ɫ��
        % ���ϴεĵ�͵�ǰ���⵽�ľ�������0.1�ף�������εĽ��������dataSizeΪǰcutLength��ԭ����Ϊ�˱���lastRange��ʼ����Ӱ��
        if  (abs(target_range(lastRangeIdx) - lastRange < 0.1) && abs(Angle_of_arrival(lastRangeIdx) - lastAoA < 5))|| dataSize < cutLength
            aa = aa + 1;
            % ����lastRange��AoA��ֵ
            lastRange = target_range(lastRangeIdx);
            lastAoA = Angle_of_arrival(lastRangeIdx);
            x_=100 * target_range(lastRangeIdx) * sind(Angle_of_arrival(lastRangeIdx));
            y_=100 * target_range(lastRangeIdx) * cosd(Angle_of_arrival(lastRangeIdx));
            
            % �������ϣ�����x��y
            mxHistory(:,dataSize) = [x_;y_];
            
            % �洢�ǶȺ;����¼
            RangeAndAoAdata(1,dataSize) = target_range(lastRangeIdx);
            RangeAndAoAdata(2,dataSize) = Angle_of_arrival(lastRangeIdx);
            if dataSize >= cutLength
                plot(mxHistory(1,dataSize), mxHistory(2,dataSize), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
                hold on;
                plot(mxHistory(1,1:dataSize), mxHistory(2,1:dataSize), 'r', 'LineWidth', 2);
                % plot(smooth(mxHistory(1,1:dataSize),31), smooth(mxHistory(2,1:dataSize),31), 'r', 'LineWidth', 2);
                set(gca,'XDir','reverse')%��X����ת
                set(gca,'YDir','reverse')%��Y����ת
                hold on;
                xlim([-limValue/2,limValue/2])
                ylim([0,limValue])
                grid on;
                xlabel('x')
                ylabel('y')
                title('tracker')
            end
        else
            bb = bb + 1;
            % ���������������εļ�¼ֵ��Ϊ�ϴμ�¼ֵ(����dataSizeһ����cutLength���ϣ���˲���Խ��)
            mxHistory(:,dataSize) = mxHistory(:,dataSize - 1);
            % �洢�ǶȺ;����¼
            RangeAndAoAdata(:,dataSize) = RangeAndAoAdata(:,dataSize - 1);
            
            plot(mxHistory(1,dataSize), mxHistory(2,dataSize), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
            hold on;
            plot(mxHistory(1,1:dataSize), mxHistory(2,1:dataSize), 'r', 'LineWidth', 2);
            % plot(smooth(mxHistory(1,1:dataSize),31), smooth(mxHistory(2,1:dataSize),31), 'r', 'LineWidth', 2);
            set(gca,'XDir','reverse')%��X����ת
            set(gca,'YDir','reverse')%��Y����ת
            hold on;
            xlim([-limValue/2, limValue/2])
            ylim([0,limValue])
            grid on;
            xlabel('x')
            ylabel('y')
            title('tracker')
        end
        %clear variables
        rvLocs=[];
        rvPks=[];
        target_signal_value=[];
        target_range=[];
        Desired_bin=[];
        iter=[];
        scPeakBin=[];
        rvPh = [] ;
        rvPhDelt = [];
        scAlphSin = [];
        Angle_of_arrival =[];
    end
    drawnow
end
save(['RangeAndAoAdata_size_',num2str(dataSize),'.mat'],'RangeAndAoAdata');
save(['mxHistory_size_',num2str(dataSize),'.mat'],'mxHistory');
save(['mxAllHistory_size_',num2str(dataSize),'.mat'],'mxAllHistory');
clearSP;