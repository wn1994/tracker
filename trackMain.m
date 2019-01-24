dbstop if error
clc
disp('******************************************************************');
addpath('..\..\RadarSystemImplementation');                           % add MATLAB API
clear all
close all                                                                      % close and delete ports
warning off
%% setup object and show properties
szPort = findRSPort;                                                              % finds port connected to board
oRS = RadarSystem(szPort);                                                        % setup object and connect to board
disp('oRS object - properties before set block:');
oRS
%% set properties
oRS.oEPRadarFMCW.lower_frequency_kHz = oRS.oEPRadarBase.min_rf_frequency_kHz+1e6; % = 57 GHz
oRS.oEPRadarFMCW.upper_frequency_kHz = oRS.oEPRadarBase.min_rf_frequency_kHz+7e6; % = 63 GHz
oRS.oEPRadarFMCW.direction = 0;                                                   % 'up-chirp'
oRS.oEPRadarFMCW.tx_power = 31;                                                   % valid powerrange: 0 - 31
% 10 will set the conducted transmit power to a value of 0 dBm
oRS.oEPRadarSoli.samplerate_Hz = 2000000;                                         % 2 MHz max. sampling rate

oRS.oEPRadarBase.num_chirps_per_frame = 31;                                        % one chirp per radar transmission
oRS.oEPRadarBase.num_samples_per_chirp = 128;                                     % 128 sample points per chirp
oRS.oEPRadarBase.rx_mask = bin2dec('1100');                                       % turn on up to 4 receiverchannels
% every bit represents one receiver, bit mask is oRS.sRXMask=[RX4 RX3 RX2 RX1]
oRS.oEPRadarSoli.tx_mode = 0;                                                     % TX1 transmitter only
% oRS.oEPRadarSoli.tx_mode = 2;                                                   % 2 for two transmitter sequencially
disp('oRS object - properties after set block:');

%% 公共参数
N=double(oRS.oEPRadarBase.num_samples_per_chirp);
Zeropad=N*4;                                                                   % = Zeropadding to multiple of N where N*4 is good choice
T=double(oRS.oEPRadarBase.chirp_duration_ns)*1e-9;                             % chirp time in seconds
BW=double(oRS.oEPRadarFMCW.upper_frequency_kHz-oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % Radar bandwidth in Hz
kf=BW/T;                                                                       % slope rate in Hz/s;
c0=3e8;                                                                        % speed of light im m/s;
disp('---------------------begin-----------------------');
%% calc_twirl Data stream初始化
global numdetections
global dopplerave
global rangeave
global angleval
global numdetectionspos
global doppleravepos
global numdetectionsneg
global doppleraveneg
BUFFER_SIZE = 100;
numdetections = zeros(1,BUFFER_SIZE);
dopplerave = zeros(1,BUFFER_SIZE);
rangeave = zeros(1,BUFFER_SIZE);
angleval = zeros(1,BUFFER_SIZE);
numdetectionspos = zeros(1,BUFFER_SIZE);
doppleravepos = zeros(1,BUFFER_SIZE);
numdetectionsneg = zeros(1,BUFFER_SIZE);
doppleraveneg = zeros(1,BUFFER_SIZE);


datavec=[];

% 上次保存的range值
lastRange = 0;
lastAoA = 0;
x=0;
y=0;
xx=0;
yy=0;
properData = 0;
improperData = 0;
% 去掉前一小段数据
cutLength = 10;

dataSize = 0;
bufferSize = 2000;
% red 线条
mxHistory = NaN(2,bufferSize);
% blue 线条
mxAllHistory = NaN(2,bufferSize);
% 存储range和AOA
RangeAndAoAdata = NaN(2,bufferSize);
hTime = figure;
set(hTime,'units','normalized','position',[0.3 0.35 0.35 0.50]);
% 坐标范围
limValue = 60;

while ishandle(hTime)
    pause(0.000001); %必须要有这个， 要不然程序可能无法得到键盘输入
    if isletter(get(hTime,'CurrentCharacter'))
        break;
    end
    
    [c, ~] = oRS.oEPRadarBase.get_frame_data;         % MATLAB API function: starting transmission on the Soli-chip and extracting raw data of every chirp, here: 40 chirps
    % 得到range-Doppler map
    [numRangeBins,numDopplerBins,doppler_FFT] = calc_range_doppler_map(c,N,Zeropad,T,BW,kf,c0,oRS);
    
    
    c = [c(:,1,1), c(:,2,1)];
    scWL=c0/(double(oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3);                    % maximum wavelength
    scSpace=0.7*scWL;                                                              % spacing between the receivers
    xscale=2*kf/c0;                                                                % scaling for the x-axis as a distance in m
    fadc =double(oRS.oEPRadarSoli.samplerate_Hz);
    xdata=linspace(0,1-1/Zeropad,Zeropad)*fadc/xscale;                             % scaling the x-axis to wproper range to frequency
    xdata1=linspace(0,1-1/N,N)*fadc/xscale;                                        % scaling for
    
    Hann_window=hann(N,'periodic');                                                % hann(Length of window, 'periodic');
    ScaleHannWin = 1/sum(Hann_window);                                             % replaces the scaling of the fft with length of spectrum with sum(window)
    
    
    % DC Removal
    % accomplished by substracting the mean of the data from the raw data suppressing the DC signal and filtering the data with a high-pass to remove all low frequencies
    % 通过从抑制DC信号的原始数据中减去数据的平均值并用高通滤波数据来消除所有低频
    
    c(:,1)=c(:,1)-mean(c(:,1));
    c(:,2)=c(:,2)-mean(c(:,2));
    Nbut=4;                                                                        % 4th order of filter
    Wn=.05;                                                                        % cutoff frequency
    [b,a]=butter(Nbut,Wn,'high');                                                  % butter(order of filter, cutoff frequency, type of filter);
    % 高通滤波
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
            %  通过取两个接收器之间的相差确定AoA，并计算AoA
            % AoA = arcsin((delta_phi/2pi)*(max. wavelength/spacing))
            
            rvPh(kk,:) = unwrap(angle(IF_info(Desired_bin(kk),:)),[]);                    % get phase information from the desired target bin for both receivers
            rvPhDelt(kk) = diff(rvPh(kk,:));                                              % get phase difference between the two receivers
            scAlphSin(kk) = (rvPhDelt(kk)/(2*pi))*(scWL/scSpace);                         % Angle of arrival formula alpha = (delta_phase/2*pi)*(lambda/0.7*lambda)
            Angle_of_arrival(kk) = asind(scAlphSin(kk));                                  % Angle value in degrees
            
            
            %% 计算神经网络需要的13个数据并得到计算结果
            detThresh=-60;
            gestureMetrics=Gesture_findNumDetections(numRangeBins,numDopplerBins,doppler_FFT,detThresh);
            % 只取检测到的最近的物体的角度
            gestureMetrics(11)=Angle_of_arrival(1);
            gestureMetrics(12)=0;
            gestureMetrics(13)=0;
            datavec=[datavec;gestureMetrics'];
            % 一个gestureMetrics大小为13，攒够4个gestureMetrics放入calc_twirl看输出结果
            if length(datavec)==13*4
                out = calc_twirl(datavec);
                datavec=[];
            end
        end
        
        %         %% plot Range angle map
        %         subplot(2,1,2)
        %         plot(target_range(1),Angle_of_arrival(1),'Marker','o','MarkerSize',10,'MarkerFaceColor','blue','MarkerEdgeColor','green');
        %         xlim([0.1,0.9])
        %         ylim([-50,50])
        %         grid on;
        %         legend(num2str(Angle_of_arrival));
        %         xlabel('Target range in m')
        %         ylabel('Azimuthal angle in degrees')
        %         title('Range angle map')
        
        %% plot tracker
        % 排除异常点
        % 找出检测到的点中离上一个点距离最近的点的index
        [~,lastRangeIdx] = min(abs(target_range - lastRange));
        
        %% 画出所有检测到的轨迹（蓝色）
%         if dataSize >= cutLength
%             xx=100 * target_range(1) * sind(Angle_of_arrival(1));
%             yy=100 * target_range(1) * cosd(Angle_of_arrival(1));
%             mxAllHistory(:,dataSize) = [xx;yy];
%             plot(xx, yy, 'or', 'MarkerFaceColor', 'blue', 'MarkerSize', 15);
%             hold on;
%             plot(mxAllHistory(1,1:dataSize), mxAllHistory(2,1:dataSize), 'blue', 'LineWidth', 2);
%             hold on;
%         end
        
        %% 画出过滤后的轨迹（红色）
        % 当上次的点和当前点检测到的距离差大于0.1米，舍弃这次的结果；放行dataSize为前cutLength的原因是为了避免lastRange初始化的影响
        if  (abs(target_range(lastRangeIdx) - lastRange < 0.04) && abs(Angle_of_arrival(lastRangeIdx) - lastAoA < 5)) || dataSize < cutLength
            properData = properData + 1;
            % 更新lastRange和AoA的值
            lastRange = target_range(lastRangeIdx);
            lastAoA = Angle_of_arrival(lastRangeIdx);
            x_=100 * target_range(lastRangeIdx) * sind(Angle_of_arrival(lastRangeIdx));
            y_=100 * target_range(lastRangeIdx) * cosd(Angle_of_arrival(lastRangeIdx));
            
            % 条件符合，更新x，y
            mxHistory(:,dataSize) = [x_;y_];
            
            % 存储角度和距离记录
            RangeAndAoAdata(1,dataSize) = target_range(lastRangeIdx);
            RangeAndAoAdata(2,dataSize) = Angle_of_arrival(lastRangeIdx);
%             if dataSize >= cutLength
%                 plot(mxHistory(1,dataSize), mxHistory(2,dataSize), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
%                 hold on;
%                 plot(mxHistory(1,1:dataSize), mxHistory(2,1:dataSize), 'r', 'LineWidth', 2);
% %                 plot(smooth(mxHistory(1,1:dataSize),31), smooth(mxHistory(2,1:dataSize),31), 'r', 'LineWidth', 2);
%                 hold on;
%                 xlim([-limValue/2,limValue/2])
%                 ylim([0,limValue])
%                 grid on;
%                 xlabel('x')
%                 ylabel('y')
%                 title('tracker')
%             end
        else
            improperData = improperData + 1;
            % 条件不符，将本次的记录值置为上次记录值(这里dataSize一定是cutLength以上，因此不会越界)
            mxHistory(:,dataSize) = mxHistory(:,dataSize - 1);
            % 存储角度和距离记录
            RangeAndAoAdata(:,dataSize) = RangeAndAoAdata(:,dataSize - 1);
            
%             plot(mxHistory(1,dataSize), mxHistory(2,dataSize), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
%             hold on;
%             plot(mxHistory(1,1:dataSize), mxHistory(2,1:dataSize), 'r', 'LineWidth', 2);
% %             plot(smooth(mxHistory(1,1:dataSize),31), smooth(mxHistory(2,1:dataSize),31), 'r', 'LineWidth', 2);
%             hold on;
%             xlim([-limValue/2, limValue/2])
%             ylim([0,limValue])
%             grid on;
%             xlabel('x')
%             ylabel('y')
%             title('tracker')
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
%     drawnow
end
clearSP;                                                                        % clear all ports