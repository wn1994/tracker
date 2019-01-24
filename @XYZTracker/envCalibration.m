function envCalibration(obj, ~, ~)

hCalFFT = figure();
hCalTime = figure();
h = waitbar(0,'Calibrating...');
szDlgDefault='Save';
mxMeanMag = zeros(obj.oRS.oEPRadarBase.num_samples_per_chirp,4);
for n=1:obj.nCalRuns
    [mxRawData, ~] = obj.oRS.oEPRadarBase.get_frame_data;

    % dispatch data as a circle starting with Rx2 up/left and going cw
    mxRawData = mxRawData-.5;
    mxFFT = fft(mxRawData);

    mxMag = abs(mxFFT);
    mxPh = angle(mxFFT);

    mxMeanMag = mxMeanMag+mxMag./obj.nCalRuns;

    mxMagCorr = mxMag - mxMeanMag;
    mxMagCorr = mxMagCorr.*(mxMagCorr>0);
    mxCorrData = real(ifft(mxMagCorr.*exp(1j*mxPh)));

    % plot time signal
    if isvalid(hCalTime)
	set(groot,'CurrentFigure',hCalTime);
	plot(mxCorrData)
    v=axis;
    axis([v(1:2) -.5 .5]);
    title('Time Data');
    end

    % plot fft signal
    if isvalid(hCalFFT)
    set(groot,'CurrentFigure',hCalFFT);
    plot(mxMagCorr)
    title('FFT Data');
    end
    
    if isvalid(h)
        waitbar(n/obj.nCalRuns)
    else
        szDlgDefault = 'Discard';
        break
    end
end

if isvalid(h)
    close(h)
end

choice = questdlg('Calibration done.', '',...
    'Discard', 'Repeat', 'Save', szDlgDefault);

if isvalid(hCalTime)
    close(hCalTime)
end
if isvalid(hCalFFT)
    close(hCalFFT)
end

switch choice
    case {'Discard' 0}
        % do nothing, as the old calibration data is still stored in the obj
    case 'Repeat'
        obj.envCalibration;
    case 'Save'
        % dispatch data as a circle starting with Rx2 up/left and going cw
        %
        % /----------\
        % | Rx2  Rx1 |
        % |          |
        % | Rx4  Rx3 |
        % +----------+
        %
        % -> Rx2, Rx1, Rx3, Rx4
        mxMeanMag = [mxMeanMag(:,2) , mxMeanMag(:,1) , mxMeanMag(:,3) , mxMeanMag(:,4)];
        save('CalibrationData.mat', 'mxMeanMag');
        obj.mxMeanMag = mxMeanMag;
end

end