classdef XYZTracker < handle

    properties
        % radar system object
        oRS;
        mxMeanMag;

        % radar data
        mxRawData;
        bAppCal = 1;
        
        % AoA algorithm properties and constants
        scNfft = 2^12; % fft length with zero padding (rows)
        scSpace = 3.5e-3; % spacing between two receivers in m (scWL*.7 for scWL=5e-3, F=6GHz)
        rvBP = [0.024 0.203]; % % normalized bandpass cutoff with 1^= nyquist
        
        % targets
        rvRange
        rvAngle
        mxTargets
        cvTargetsMean
        nTaps = 40;
        mxHistory;
        
        % gestures
        cvGesture = NaN(3,1);
        
        % calibration
        nCalRuns = 50;
        
        % 过滤无效数据
        lastRange = 0;
        i = 0;
    end

    methods
        % constructor
        function obj = XYZTracker(szPort)
            % setup radar system object oRS
            if ~nargin % hand over known COM port in szPort (e.g. 'COM3') to speed up setup
                szPort = findRSPort; % scan all available ports
            end

            obj.oRS = RadarSystem(szPort); % setup object and connect to board

            if isprop(obj.oRS, 'oEPRadarSoliC') % check for SoliC
                obj.oRS.oEPRadarSoliC.repeat_chip_setup
                obj.oRS.oEPRadarSoliC.enable_easy_mode(1);
            end

            % set parameters
            obj.oRS.oEPRadarFMCW.lower_frequency_kHz = 57000000; % kHz
            obj.oRS.oEPRadarFMCW.upper_frequency_kHz = 63500000; % kHz
            obj.oRS.oEPRadarFMCW.direction = 'Up Only';
            obj.oRS.oEPRadarFMCW.tx_power = obj.oRS.oEPRadarBase.max_tx_power;
            obj.oRS.oEPRadarSoli.samplerate_Hz = 3000000;
            obj.oRS.oEPRadarBase.num_chirps_per_frame = 1;
            obj.oRS.oEPRadarBase.num_samples_per_chirp = 256;
            obj.oRS.oEPRadarBase.rx_mask = bin2dec('1111');
            obj.oRS.oEPRadarSoli.tx_mode = 'Only Tx1';
            disp('Radar System object')
            disp(obj.oRS)

            % load free field filters or local calibration data if available
            obj.mxMeanMag = zeros(obj.oRS.oEPRadarBase.num_samples_per_chirp,4);
            
            obj.lastRange = 0;
            obj.i = 0;
            try % either general free field data...
                s = load('RawFilter.mat', 'mxMeanMag');
                obj.mxMeanMag = s.mxMeanMag;
                clear s;
            catch
            end

            try % ...or preferable calibration data over free field filters
                s = load('CalibrationData.mat', 'mxMeanMag');
                obj.mxMeanMag = s.mxMeanMag;
                clear s;
            catch
            end
            
            % setup history
            obj.mxHistory = NaN(3,obj.nTaps);

        end

        function getFrameData(obj)
            % get data from chip
            [mxRawData_, ~] = obj.oRS.oEPRadarBase.get_frame_data;
            mxRawData_ = mxRawData_-.5; % subtract offset

            % dispatch data as a circle starting with Rx2 up/left and going clockwiese
            %
            % /----------\
            % | Rx2  1  Rx1 |
            % |   4       |    2
            % | Rx4  3  Rx3 |
            % +----------+
            %
            % -> Rx2, Rx1, Rx3, Rx4
            obj.mxRawData = [mxRawData_(:,2) , mxRawData_(:,1) , mxRawData_(:,3) , mxRawData_(:,4)];
            % apply calibration data if desired
            if obj.bAppCal
                % build FFT
                mxFFT = fft(obj.mxRawData);
                mxMag = abs(mxFFT);
                mxPh = angle(mxFFT);

                % subtract stored mean magnitude and rectify
                mxMagCorr = mxMag - obj.mxMeanMag;
                mxMagCorr = mxMagCorr.*(mxMagCorr>0);
%                 real()取复数的实部
%                 mxRawData  维度256x4
                obj.mxRawData = real(ifft(mxMagCorr.*exp(1j*mxPh)));
            end            
        end
        
        function findTargets(obj)
            % properties and constants
            scF0 = double(obj.oRS.oEPRadarFMCW.lower_frequency_kHz)*1e3; % start freq in Hz (57e9)
            scF1 = double(obj.oRS.oEPRadarFMCW.upper_frequency_kHz)*1e3; % start freq in Hz (63.5e9)            
            scTsw = double(obj.oRS.oEPRadarBase.chirp_duration_ns)*1e-9; % sweep time in s (150e-6)
            

            % get angle of arrival and range of the target
            [obj.rvRange, obj.rvAngle, ~] = AoADetection(obj.mxRawData, scF0, scF1, scTsw, obj.scSpace, obj.scNfft, obj.rvBP);
             obj.i =  obj.i + 1;
             if  abs(obj.rvRange(1) - obj.lastRange < 0.05) || obj.i < 50
                 obj.lastRange = obj.rvRange(1);
            % For a uniform circular array the positive angle values refer
            % to a counter clockwiese direction. Invert values 1 and 4 to
            % grow positive from bottom to top and left to right.
            % 对于一个均匀的圆形阵列，正的角度值指的是逆时钟方向。颠倒数值1和4，使其从下到上、从左到右逐渐变为正数。
            %  before        :       after
            %  A <-- B       :       A --> B
            %  |     ^       :       ^     ^
            %  v     |       :       |     |
            %  D --> C       :       D --> C           
            obj.rvAngle([1 4]) = -obj.rvAngle([1 4]);

            % build xyz coordinates for each Rx
%            AoA2xyz 得到每个Rx的位置xyz
            obj.mxTargets = AoA2xyz(obj.rvRange, obj.rvAngle);

            % build averages
%             cvTargetsMean   3x1
            obj.cvTargetsMean = mean(obj.mxTargets,2);
            disp(obj.mxTargets(:,1));
            % check for mirrored/erroneos targets
%             cvTargetsStd = std(obj.mxTargets,1,2);
%             scStdSum = sum(cvTargetsStd(1:2).^2);            
%             if (scStdSum>0.025)
%                 obj.cvTargetsMean = NaN(3,1);
%             end

            % update history
%             mxHistory 3x40
% mxHistory列每次右移一位，在第一位存放当次的obj.cvTargetsMean
            obj.mxHistory(:,2:end) = obj.mxHistory(:,1:end-1); % shift right
            obj.mxHistory(:,1) = obj.mxTargets(:,1);;
             end
        end

        function detectGestures(obj)
            % compute mean velocity
            scVTap = 5;
%             求路径的一阶导数，得到速度
% diff(obj.mxHistory(:,1:scVTap),1,2)，计算obj.mxHistory(:,1:scVTap)列之间的差分（第n列减去第n-1列）
% mean(X,2)计算每行的均值
            cvVelocity = -mean(diff(obj.mxHistory(:,1:scVTap),1,2),2);
            scXYGestThres = 0.02; % m
            scZGestThres = 0.02; % m

            % xy gestures
%             cart2pol直角坐标转变为极坐标
%             scTheta是沿逆时针方向与X轴正方向的夹角，scR投影与原点的距离 
            [scTheta, scR] = cart2pol(cvVelocity(1), cvVelocity(2));
%             划定速度范围,
            if (scR>scXYGestThres)&&(scR<.1) % (scR<.1) to avoid gesture triggers from errors (flips)
%                 将方向固定在了45度或者坐标轴的方向
                scTheta = .25*pi*round(scTheta/(.25*pi));
                [obj.cvGesture(1), obj.cvGesture(2)] = pol2cart(scTheta, scR);
            else
                obj.cvGesture(1:2) = NaN(2,1);
            end
%  obj.cvGesture(1:2)
            % z gestures
            if (abs(cvVelocity(3))>scZGestThres)&&(scR<.2) % (scR<.2) to avoid gesture triggers from errors (flips)
                obj.cvGesture(3) = cvVelocity(3);
            else
                obj.cvGesture(3) = NaN;
            end
        end

        % envelope calibration method
        envCalibration(obj, ~, ~)
    end

end
