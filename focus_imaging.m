%% 相控阵探头参数
f0 = 4.0e6;              % 中心频率 3.5 MHz
element_num = 64;        % 阵元数量
height = 5e-3;           % 阵元高度 (m)
width = 0.12e-3;         % 阵元宽度 (m)
kerf = 0.18e-3;          % 阵元间距 (m)
pitch = width + kerf;
focus = 60e-3;           % 发射聚焦深度 (m)
c = 1540;                % 声速 (m/s)
fs = 25e6;              % 采样频率 (Hz)
fc = 1.5e6;           % 滤波截至凭（Hz）
wvln = c / f0;      % 波长
Ts = 1/fs;               % 采样间隔 (s)

% 阵元参数
x_ele = ([0:element_num-1]-(element_num-1)/2).*pitch;
z_ele  = zeros(1,length(x_ele));

% Imaging Parameters
num_MLA = 1;
num_Scan = 64;
fov_Scan = 90;
depth = 4096;
F = 85e-3;


%% 扇扫参数设置
% scan_angle_deg = -30:2:30;
angle_tx = fov_Scan - fov_Scan/(num_Scan*num_MLA)*((num_MLA/2-1)*2+1);
angles = linspace(-angle_tx/2, angle_tx/2, num_Scan);
theta = deg2rad(angles);
num_line = length(angles);
raw_data = cell(1, num_line); % 存储每条扫描线的多阵元原始数据
tstart = zeros(num_line,1);


for i = 1:64
    rfdata = readmatrix(['rfdata\rfdata_1_', num2str(i), '.csv']);
%     rfdata = (rfdata - 512) ;
    
    rfdata = bandpass_filter(rfdata);
%     rfdata = rfdata - mean(rfdata(:));
    
    raw_data{i} = rfdata;
    tstart(i) = 0;
end

%% 接收波束合成

rx_num_line = num_Scan*num_MLA;
ele_x = (0:element_num-1)*pitch-(element_num-1)/2.*pitch;
ele_y = zeros(size(ele_x));
ele_z = zeros(size(ele_x));
ele_xyz = [ele_x', ele_y', ele_z'];

%% 参数预定义
das = zeros(rx_num_line, depth);
grid_x = zeros(rx_num_line, depth);
grid_y = zeros(rx_num_line, depth);

% 变迹
apod_win = Beamform.dynamic_apodization(element_num,depth,c,fs, pitch, 0.8, [6,64], 'chebwin',50);

%% 波束合成浮点计算
step_angle = angles(2) - angles(1);
% step_list = [-0.375 -]
for i = 1 : num_Scan
    for k = 1:num_MLA
        angle = angles(i) + step_angle * (k - (num_MLA+1)/2) / num_MLA;
        rfdata = raw_data{i}';
        
        % 延时计算
        Fj_float = (0:1:depth-1) * (0.5*c/fs);                  % 每个深度点实际距离
        line_x = ele_x' - Fj_float .* sind(angle);              % 每个阵元与焦点的横向差值
        line_y = ele_y' - Fj_float .* cosd(angle);              % 每个阵元与焦点的纵向差值
        dist_float = sqrt(line_x.^2 + line_y.^2) + Fj_float;    % 发射+接收距离         
        
        %     delay_float = dist_float ./ (c/fs);
        delay_float = (dist_float / c - tstart(i)) * fs;

        % DAS (delay and sum)
        xc = 1 : size(rfdata, 2);
        for j = 1 : element_num
            foc(j, :) = interp1(xc, rfdata(j, :), delay_float(j, :), 'linear', 0.0);
        end
        
        
        foc = foc .* apod_win;
        beamdata = sum(foc);
        
        % IQ解调滤波
        m = size(beamdata,2);
        w = 2 * pi * f0;
        x = (0:m-1) / fs;
        I = cos(w*x).*beamdata;
        Q = sin(w*x).*beamdata;

        % fir滤波
        n = 63;
        wn = fc / (fs/2); % 归一化截至频率
        fir = fir1(n, wn, 'low');
        
        If = filter(fir, 1, I);
        Qf = filter(fir, 1, Q);
        
        % 包络检波
        Fdata = sqrt(If.^2 + Qf.^2);
        
        das((i-1)*num_MLA + k, :) = Fdata;

        grid_x((i-1)*num_MLA + k, :) = Fj_float .* sind(angle);
        grid_y((i-1)*num_MLA + k, :) = Fj_float .* cosd(angle);
    end
end


%% scan convert
rmax = 158e-3;
rlims = [0, rmax];
xlims = rlims(2) * [-0.7, 0.7];
zlims = rlims(2) * [0, 1];

% sector_half_angle = 45;  % 半开角45度
% xlims = rmax * [-sind(sector_half_angle), sind(sector_half_angle)];
% make pixel grid (Cartesian pixel grid)
x = xlims(1):wvln:xlims(2);
z = zlims(1):wvln:zlims(2);

[xx, zz] = meshgrid(x, z);
yy = zeros(size(xx));

img_grid = cat(3, xx, yy, zz);


grid_x = grid_x(:);
grid_y = grid_y(:);
img_grid_x = img_grid(:, :, 1);
img_grid_x = img_grid_x(:);
img_grid_y = img_grid(:, :, 3);
img_grid_y = img_grid_y(:);

bimgsc = griddata(grid_x, grid_y, das(:), img_grid_x, img_grid_y, 'linear');
bimgsc(isnan(bimgsc)) = 1e-22;
bimg = reshape(bimgsc, size(img_grid, 1), size(img_grid, 2));

drange = 60;
bimg = abs(bimg);
bimg = 20 * log10(bimg);
bimg = bimg - max(bimg(:));
extent = [img_grid(1, 1, 1), img_grid(1, end, 1), img_grid(1, 1, 3), img_grid(end, 1, 3)] * 1e3;
gcf = figure;
imagesc(extent([1 2]), extent([3 4]), bimg);
colormap gray;
caxis([-drange 0]);
% set(gca, 'YDir', 'normal');
xlabel('Lateral distance [mm]');
ylabel('Axis distance [mm]');
axis image;
print(gcf, 'image_high_quality.jpg', '-djpeg', '-r300');





function y = bandpass_filter(x)

    % Fs = 20e6;             % 采样频率 20MHz
    % Fstop1 = 0.9e6;        % 下阻带截止频率 0.9MHz
    % Fpass1 = 1e6;          % 通带下限频率 1MHz
    % Fpass2 = 5e6;          % 通带上限频率 5MHz
    % Fstop2 = 5.1e6;        % 上阻带截止频率 5.1MHz
    % Astop1 = 60;           % 下阻带衰减 60dB
    % Apass = 1;             % 通带波纹 1dB
    % Astop2 = 60;           % 上阻带衰减 60dB

    % d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',...
    % Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
    % Hd = design(d, 'equiripple');

    % y = filter(Hd, x);

    Fs = 25e6;  % 20 MHz

    fpass = [1e6 10e6];  % 1 MHz - 5 MHz
    wn = fpass / (Fs / 2);

    N = 128;  % 可以调整阶数
    b = fir1(N, wn, 'bandpass', hamming(N+1));
%     b = b/sum(b(:));
    y = zeros(size(x));
    tap_half = N/2;
    
    for i = 1:size(x,2)
         temp = filter(b, 1, [zeros(tap_half,1); x(:,i); zeros(tap_half,1)]);
         y(:,i) = temp(N+1:end,1);
    end
%     figure;plot( x(:,1));hold on;plot(y(:,1));legend('Oroginal','filterd')
end

