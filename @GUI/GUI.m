classdef GUI < handle

    properties
        hGUI
        handles
        hTime
        hFFT
        oXYZ
    end

    methods
        % constructor
        function obj = GUI(oXYZ)
            % store handle to obj
            obj.oXYZ = oXYZ;

            % load and init .fig
            obj.hGUI = openfig('GUI','new','invisible');
            obj.handles = guihandles(obj.hGUI);

            % bind callbacks
            obj.handles.menuEnv.Callback = @oXYZ.envCalibration;
            obj.handles.menuAppCal.Callback = @obj.applyCal;
            obj.handles.menuGestures.Callback = @obj.showGestures;
            obj.handles.menuRawData.Callback = @obj.showRawData;

            % load logo and set alpha channel
            [A, map, alpha] = imread('Infineon_logo.png');
            h = imshow(A, map, 'Parent', obj.handles.axesLogo);
            set(h, 'AlphaData', alpha);

            % clear axes and draw xy axis
            obj.clearAxes;

            % switch on initialized GUI
            set(obj.hGUI,'Visible', 'on');
        end

        % "Apply Calibration Data" "Ctrl+A"
        function applyCal(obj, hMenu, ~)
            obj.oXYZ.bAppCal = ~obj.oXYZ.bAppCal;
            if obj.oXYZ.bAppCal
                hMenu.Checked = 'on';
            else
                hMenu.Checked = 'off';
            end
        end

        % "Show Gestures" "Ctrl+G" callback
        function showGestures(~, hMenu, ~)
            if strcmp(hMenu.Checked, 'off')
                hMenu.Checked = 'on';
            else
                hMenu.Checked = 'off';
            end
        end

        % "Show Raw Data" "Ctrl+R" callback
        function showRawData(obj, hMenu, ~)
            if strcmp(hMenu.Checked, 'off')
                hMenu.Checked = 'on';
                obj.hTime = figure;
                obj.hFFT = figure;                
            else
                hMenu.Checked = 'off';
                if isvalid(obj.hTime)
                    close(obj.hTime);
                end
                if isvalid(obj.hFFT)
                    close(obj.hFFT);
                end
            end
        end

        % clear axes
        function clearAxes(obj)
            % xy axes
            cla(obj.handles.axesXY)
            hold(obj.handles.axesXY,'on');
            v = axis(obj.handles.axesXY);
            line(v(1:2), [0 0], 'Color', 'k', 'Parent', obj.handles.axesXY);
            line([0 0], v(3:4), 'Color', 'k', 'Parent', obj.handles.axesXY);
            
            % z axis
            cla(obj.handles.axesZ);
            axis(obj.handles.axesZ, [.45 0.55 0 0.6])
        end

        % plot xyz targets
        function plotTargets(obj)
            set(groot,'CurrentFigure',obj.hGUI);

            set(obj.hGUI,'CurrentAxes', obj.handles.axesXY)
            plot(obj.oXYZ.cvTargetsMean(1), obj.oXYZ.cvTargetsMean(2), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 12)
            plot(obj.oXYZ.mxHistory(1,:), obj.oXYZ.mxHistory(2,:),'r','LineWidth', 2);

            % plot heigths
            set(obj.hGUI,'CurrentAxes', obj.handles.axesZ)
            line([0.455 .545],[obj.oXYZ.cvTargetsMean(3) obj.oXYZ.cvTargetsMean(3)], 'Color', 'r', 'LineWidth', 5);

            % disp output
            if isnan(obj.oXYZ.cvTargetsMean)
                obj.handles.textTarget.String = '';
            else
                obj.handles.textTarget.String = sprintf('Target (x y z): %1.2fcm    %1.2fcm    %1.2fcm', obj.oXYZ.cvTargetsMean*100);
            end
        end

        % plot raw data
        function plotRawData(obj)
            set(groot,'CurrentFigure',obj.hGUI);
            set(obj.hGUI,'CurrentAxes',obj.handles.axesXY)
            
            % plot x/y
            set(obj.hGUI,'CurrentAxes',obj.handles.axesXY)
            
            text(obj.oXYZ.mxTargets(1,1), obj.oXYZ.mxTargets(2,1), 'RX2', 'Color', [.5 .5 .5 .5], 'HorizontalAlignment', 'center', 'FontSize', 7)
            text(obj.oXYZ.mxTargets(1,2), obj.oXYZ.mxTargets(2,2), 'RX1', 'Color', [.5 .5 .5 .5], 'HorizontalAlignment', 'center', 'FontSize', 7)
            text(obj.oXYZ.mxTargets(1,3), obj.oXYZ.mxTargets(2,3), 'RX3', 'Color', [.5 .5 .5 .5], 'HorizontalAlignment', 'center', 'FontSize', 7)
            text(obj.oXYZ.mxTargets(1,4), obj.oXYZ.mxTargets(2,4), 'RX4', 'Color', [.5 .5 .5 .5], 'HorizontalAlignment', 'center', 'FontSize', 7)
            
            % plot heigths
            set(obj.hGUI,'CurrentAxes',obj.handles.axesZ)
            
            line([0 1],[obj.oXYZ.mxTargets(3,1) obj.oXYZ.mxTargets(3,1)], 'Color', [.5 .5 .5 .5]);
            line([0 1],[obj.oXYZ.mxTargets(3,2) obj.oXYZ.mxTargets(3,2)], 'Color', [.5 .5 .5 .5]);
            line([0 1],[obj.oXYZ.mxTargets(3,3) obj.oXYZ.mxTargets(3,3)], 'Color', [.5 .5 .5 .5]);
            line([0 1],[obj.oXYZ.mxTargets(3,4) obj.oXYZ.mxTargets(3,4)], 'Color', [.5 .5 .5 .5]);
            
            % output strings
            sRange=sprintf('sRange\n   %2.1f cm   %2.1f cm\n%2.1f cm    %2.1f cm', obj.oXYZ.rvRange([1 2 4 3])*100);
            sAngle=sprintf('sAngle\n   %2.1f \n %2.1f    %2.1f \n	%2.1f', obj.oXYZ.rvAngle([1 4 2 3])*360./(2*pi)); % convert from rad to deg

            % plot the data
            if isvalid(obj.hTime)
                set(groot,'CurrentFigure',obj.hTime);
                plot(real(obj.oXYZ.mxRawData))
                v=axis;
                axis([v(1:2) -.5 .5]);
                title(sprintf('%s\n%s',sRange,sAngle));
            end
            
            % plot the data
            if isvalid(obj.hFFT)
                set(groot,'CurrentFigure',obj.hFFT);
                plot(abs(fft(obj.oXYZ.mxRawData(1:round(end/2),:))));
                title(sprintf('%s\n%s',sRange,sAngle));
            end            
        end

        % plot gestures
        function plotGestures(obj)
            set(groot,'CurrentFigure',obj.hGUI);
            
            sf = 5; % scale factor
            % plot x/y
            if(~isnan(obj.oXYZ.cvGesture(1:2)))
                set(obj.hGUI,'CurrentAxes',obj.handles.axesXY)
                quiver(-sf*obj.oXYZ.cvGesture(1), -sf*obj.oXYZ.cvGesture(2), sf*2*obj.oXYZ.cvGesture(1), sf*2*obj.oXYZ.cvGesture(2), 1, 'LineWidth', 8, 'MaxHeadSize', .5 ,'Color', 'b')
            end
            
            % plot heigths
            if(~isnan(obj.oXYZ.cvGesture(3)))
                set(obj.hGUI,'CurrentAxes',obj.handles.axesZ)
                quiver(.5, .3-sf*obj.oXYZ.cvGesture(3), 0, sf*2*obj.oXYZ.cvGesture(3), 1, 'LineWidth', 8, 'MaxHeadSize', .5 ,'Color', 'b')
            end
        end
    end

end
