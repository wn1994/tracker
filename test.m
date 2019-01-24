figure;
i = 0;
while true %遗传算法的大循环部分
    % 遗传算法优化
    
    % 更新图形, 用于查看优化效果
    
    % 点击图形界面， 然后点击任意一个字母按键
    pause(0.000001); %必须要有这个， 要不然程序可能无法得到你的键盘输入
    i = i + 1
    if isletter(get(gcf,'CurrentCharacter'))
        break;
    end
end

% 继续做其他事
disp('继续做了其他事');