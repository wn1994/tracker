% data
dataSize = 185;
RangeAndAoAdata = struct2array(load(['RangeAndAoAdata_size_',num2str(dataSize),'.mat']));
target_range = RangeAndAoAdata(1,:);
Angle_of_arrival = RangeAndAoAdata(2,:);

% 去掉前一小段数据
cutLength = 15;
% 坐标范围
limValue = 60;
% 是否做平滑处理
isSmooth = false;

if isSmooth
    target_range = smooth(target_range,31);
    Angle_of_arrival = smooth(Angle_of_arrival,31);
end

% 原始数据
trackOrigin = NaN(2,bufferSize);
% 对AoA进行/2平均的数据
track_2 = NaN(2,bufferSize);
% 对AoA进行/3平均的数据
track_3 = NaN(2,bufferSize);

for i = 3:dataSize
    x_=100 * target_range(i) * sind(Angle_of_arrival(i)) - 10;
    y_=100 * target_range(i) * cosd(Angle_of_arrival(i));
    trackOrigin(:,i) = [x_;y_];
    calc_AoA = (Angle_of_arrival(i)+Angle_of_arrival(i - 1))/2;
    x_=100 * target_range(i) * sind(calc_AoA) - 5;
    y_=100 * target_range(i) * cosd(calc_AoA);
    track_2(:,i) = [x_;y_];
    calc_AoA = (Angle_of_arrival(i) + Angle_of_arrival(i - 1) + Angle_of_arrival(i - 2)) / 3;
    x_=100 * target_range(i) * sind(calc_AoA);
    y_=100 * target_range(i) * cosd(calc_AoA);
    track_3(:,i) = [x_;y_];
end
% track(1,:) = smooth(track(1,:),29);
plot(trackOrigin(1,cutLength:dataSize), trackOrigin(2,cutLength:dataSize), 'red', 'LineWidth', 2);
hold on;
plot(track_2(1,cutLength:dataSize), track_2(2,cutLength:dataSize), 'blue', 'LineWidth', 2);
hold on;
plot(track_3(1,cutLength:dataSize), track_3(2,cutLength:dataSize), 'black', 'LineWidth', 2);
hold on;
xlim([-limValue/2,limValue/2])
ylim([0,limValue])
grid on;
xlabel('x')
ylabel('y')
title('tracker')
% legend('原始数据','修正数据/2','修正数据/3')