% This simple example demos the usage of the Matlab Sensing Interface
% It currently runs with the BGT60TR24 SoliB and SoliC board.

%% cleanup and init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
disp('******************************************************************');
addpath('..\..\RadarSystemImplementation'); % add Matlab API
clear all %#ok<CLSCR>
close all
resetRS; % close and delete ports


%% setup GUI and XYZTracker main object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oXYZ = XYZTracker;
oGUI = GUI(oXYZ);


%% process loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ishandle(oGUI.hGUI)

    % get and dispatch data from chip
    oXYZ.getFrameData;

    % find targets using an fft based angle of arrival and distance estimation
    oXYZ.findTargets;

    % clear main axes
    oGUI.clearAxes;
%     根据四个Rx的均值确定位置
    % plot xyz target data
    oGUI.plotTargets;

    % plot raw data
    if strcmp(oGUI.handles.menuRawData.Checked, 'on')
        oGUI.plotRawData;
    end

    % compute and plot gestures
    if strcmp(oGUI.handles.menuGestures.Checked, 'on')
        oXYZ.detectGestures;
        oGUI.plotGestures;
    end

    drawnow
end % end loop



%% cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all %#ok<CLSCR>
