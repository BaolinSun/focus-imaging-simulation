% Field II Phase Array B-mode
clear all; close all;

%% Init Field II
path(path, 'D:\MyProjects\matlab\Field_II_ver_3_30_windows');
field_init(-1);

%% 相控阵探头参数
probe = Probe();
probe.type = 'phase aray';
f0 = 2.5e6;              % 中心频率 3.5 MHz
element_num = 64;        % 阵元数量
height = 5e-3;           % 阵元高度 (m)
width = 0.12e-3;         % 阵元宽度 (m)
kerf = 0.18e-3;          % 阵元间距 (m)
pitch = width + kerf;
focus = 60e-3;           % 发射聚焦深度 (m)
c = 1540;                % 声速 (m/s)
fs = 20e6;              % 采样频率 (Hz)
Ts = 1/fs;               % 采样间隔 (s)

% 阵元参数
x_ele = ([0:element_num-1]-(element_num-1)/2).*pitch;
z_ele  = zeros(1,length(x_ele));
probe.N_elements = element_num;
probe.pitch = pitch;
probe.x_ele = x_ele;
probe.y_ele = z_ele;
probe.z_ele = z_ele;
probe.ele_pos = [x_ele; z_ele; z_ele]';

%% 创建相控阵探头（发射和接收分开定义）
% 发射探头：设置发射聚焦和偏转
Th = xdc_linear_array(element_num, width, height, kerf, 1, 10, [0 0 focus]);
% 接收探头：禁用接收聚焦，独立阵元接收
Rh = xdc_linear_array(element_num, width, height, kerf, 1, 10, [0 0 0]);


% 设置2个周期高斯脉冲相应、1个周期激励脉冲
dt  = 1/fs;
t0 = (-1/f0): dt:(1/f0);
impulse_response = gauspuls(t0, f0);
impulse_response = impulse_response-mean(impulse_response);
pulse_duration = 1;
te = 0:dt:pulse_duration/f0;
excitation = square(2*pi*f0*te);
% 设置激励脉冲
xdc_excitation (Th, excitation);
% 设置脉冲相应
xdc_impulse (Th, impulse_response);
xdc_impulse (Rh, impulse_response);

%% 设置全局参数
set_sampling(fs);
set_field('c', c);

%% 生成散射点
point_position(1,:) = [0 0 20e-3];
point_position(2,:) = [0 0 30e-3];
point_position(3,:) = [-5e-3 0 30e-3];
point_position(4,:) = [5e-3 0 30e-3];
point_position(5,:) = [0 0 40e-3];
point_position(5,:) = [0 0 150e-3];
point_amplitudes = ones(size(point_position,1),1);

%% 扇扫参数设置
F = 120e-3;
% scan_angle_deg = -30:2:30;
scan_angle_deg = linspace(-45, 45, 64);
num_lines = length(scan_angle_deg);
rf_data_multi = cell(1, num_lines); % 存储每条扫描线的多阵元原始数据

%% 主循环：逐角度发射，记录各阵元原始回波
for line = 1:num_lines
    % --- 发射设置：聚焦与偏转 ---
    angle_deg = scan_angle_deg(line);
    angle_rad = angle_deg * pi/180;
    
    % 发射延时
    emit_delay = phase_array_transmit_delay(probe, angle_rad, F, c);
    
    xdc_apodization(Th, 0, ones(1, element_num)); % 发射孔径全开
    xdc_focus_times(Th, 0, emit_delay); % 应用发射延迟
    
    % --- 接收设置：禁用聚焦，获取所有阵元原始信号 ---
    xdc_apodization(Rh, 0, ones(1, element_num));  % 接收孔径全开
    xdc_focus_times(Rh, 0, zeros(1, element_num)); % 接收延迟设为0（禁用聚焦）
    
    [rf_multi, t_start] = calc_scat_multi(Th, Rh, point_position, point_amplitudes);
    % rf_multi 结构：[时间采样点 × 接收阵元]
    
    % 存储数据（需确保时间轴对齐）
    rf_data_multi{line} = rf_multi;
end

%% 关闭探头并清理内存
xdc_free(Th);
xdc_free(Rh);
field_end;

%% 保存rfdata为csv文件
for i = 1:num_lines
    rfdata = rf_data_multi{i};
    writematrix(rfdata, ['rfdata\rfdata_1_', num2str(i), '.csv'])
end


%% 接收波束合成






%% 自定义局部函数

function grid = make_foctx_grid(rlims, dr, dirs)

    r = rlims(1) : dr : rlims(2);
    t = dirs(:, 1);
    [tt, rr] = meshgrid(t, r);

    xx = rr .* sin(tt);
    zz = rr .* cos(tt);
    yy = zeros(size(xx));

    grid = cat(3, xx, yy, zz);
end
