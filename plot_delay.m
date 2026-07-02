N = 64;
Fs = 100e6;
c = 1540;
pitch = 0.254e-3;

fc_depth = 200e-3;

num_line = 64;
angles = linspace(-45, 45, 64);
theta = deg2rad(angles);

figure;
for line = 1:num_line
    angle_rad = theta(line);    
    focus_point = [angle_rad,fc_depth];

    [delay_s,delay_clk] = us_calc_tx_delay(N,Fs,c,pitch,focus_point);
    
    plot(delay_s);
    hold on;
end


addpath("probe")

%% 相控阵探头参数
probe = Probe('phase array');
f0 = 3e6;              % 中心频率 3.5 MHz
c = 1540;
element_num = 64;        % 阵元数量
pitch = 0.254e-3;
focus = 60e-3;           % 发射聚焦深度 (m)


% 阵元参数
x_ele = ([0:element_num-1]-(element_num-1)/2).*pitch;
z_ele  = zeros(1,length(x_ele));
probe.N_elements = element_num;
probe.pitch = pitch;
probe.x_ele = x_ele;
probe.y_ele = z_ele;
probe.z_ele = z_ele;
probe.ele_pos = [x_ele; z_ele; z_ele]';

F = 200e-3;
angles = linspace(-45, 45, 64);
theta = deg2rad(angles);
figure;
for line = 1:num_line
    % --- 发射设置：聚焦与偏转 ---
    angle_rad = theta(line);
    
    % 发射延时
    emit_delay = phase_array_transmit_delay(probe, angle_rad, F, c);
    plot(emit_delay);
    hold on;
end


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


function emit_delay = phase_array_transmit_delay(probe, angle_rad, F, c)
    
    xx = F * sin(angle_rad);
    zz = F * cos(angle_rad);
    yy = 0 * xx;

    focus_point = [xx; yy; zz];
    focus_point = focus_point';

    aperture_positions = [probe.x_ele; probe.y_ele; probe.z_ele];
    aperture_positions = aperture_positions';

    % Calculate distances from each element to the focus point
    L = sqrt(sum((aperture_positions - focus_point).^2, 2));

    emit_delay = 5e-6 - (L - F) / c;
    emit_delay = emit_delay';
end