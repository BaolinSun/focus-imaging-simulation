function [delay_s,delay_clk] = us_calc_tx_delay(N,Fs,c,pitch,focus_point)
    %US_CALC_TX_DELAY_THETAR  TX delay for a linear array to focus/steer to (theta, R)
    % focus = [theta(rad), R(m)]
    %
    % 输出 delay_s 的定义：每个阵元相对最早发射阵元的额外延时（最小值为 0）
    % 这与 Field II 的 xdc_focus_times(Th,0,delay_s) 的用法最自洽。
    
    x  = 0:N-1;
    xi = (x - (N-1)/2) * pitch;      % 阵元中心坐标（以阵列中心为0）
    theta = focus_point(1);
    fc_depth = focus_point(2);
    
    x_f = fc_depth*sin(theta);
    z_f = fc_depth*cos(theta);
    
    t = sqrt((xi - x_f).^2 + z_f.^2) / c;   % 各阵元到焦点的传播时间
    
    delay_s = fc_depth/c - t;
    % 关键：让最早发射的阵元延时为0，避免任何“魔法偏置”
    delay_s = delay_s - min(delay_s);
    
    
    delay_clk = round(delay_s * Fs);
end
