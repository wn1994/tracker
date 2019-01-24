% data
dataSize = 185;
mxHistory = struct2array(load(['mxHistory_size_',num2str(dataSize),'.mat']));
mxAllHistory = struct2array(load(['mxAllHistory_size_',num2str(dataSize),'.mat']));
RangeAndAoAdata = struct2array(load(['RangeAndAoAdata_size_',num2str(dataSize),'.mat']));
% 去掉前一小段数据
cutLength = 15;
% 拖尾历史数据长度
tail = 20;
% 坐标范围
limValue = 60;
% for i = tail + 1:dataSize
%     if i ==  tail + 1
%         pause(1);
%     end
%     plot(mxHistory(1,i), mxHistory(2,i), 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
%     hold on;\       
%     plot(mxHistory(1,i:-1:i - tail), mxHistory(2,i:-1:i - tail),'r','LineWidth', 2);
%     hold on;
%     plot(mxAllHistory(1,i), mxAllHistory(2,i), 'or', 'MarkerFaceColor', 'blue', 'MarkerSize', 12);
%     hold on;
%     plot(mxAllHistory(1,i:-1:i - tail), mxAllHistory(1,i:-1:i - tail),'blue','LineWidth', 2);
%     hold on;
%     xlim([-limValue/2,limValue/2])
%     ylim([0,limValue])
%     grid on;
%     xlabel('x')
%     ylabel('y')
%     title('tracker')
%     drawnow
%     clf
% end
pause(0.2);
subplot(2,2,1)
plot(mxAllHistory(1,cutLength:dataSize), mxAllHistory(2,cutLength:dataSize), 'blue', 'LineWidth', 2);
set(gca,'XDir','reverse')%对X方向反转
set(gca,'YDir','reverse')%对Y方向反转
hold on;
xlim([-limValue/2,limValue/2])
ylim([0,limValue])
grid on;
xlabel('x')
ylabel('y')
title('tracker')

subplot(2,2,2)
plot(mxHistory(1,cutLength:dataSize), mxHistory(2,cutLength:dataSize), 'red', 'LineWidth', 2);
set(gca,'XDir','reverse')%对X方向反转
set(gca,'YDir','reverse')%对Y方向反转
hold on;
xlim([-limValue/2,limValue/2])
ylim([0,limValue])
grid on;
xlabel('x')
ylabel('y')
title('tracker')

subplot(2,2,3)
plot([cutLength:dataSize], RangeAndAoAdata(1,cutLength:dataSize), 'black', 'LineWidth', 1);
hold on;
grid on;
xlabel('frame')
ylabel('range')
title('range')

subplot(2,2,4)
plot([cutLength:dataSize], RangeAndAoAdata(2,cutLength:dataSize), 'black', 'LineWidth', 1);
hold on;
grid on;
xlabel('frame')
ylabel('AoA')
title('AoA')



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