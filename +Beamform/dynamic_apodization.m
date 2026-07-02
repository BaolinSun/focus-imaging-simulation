function [apod_weights, active_info] = dynamic_apodization(nElement,nSample,c,fs, pitch, rx_fnumber, min_max_aper, win_type, varargin)
% DYNAMIC_APODIZATION 生成动态变迹权重表
%
% 输入:
%   channel_data  - 通道数据 [nElement, nSample]，行是阵元，列是采样点
%   rx_fnumber    - 接收f数，用于控制深度与孔径的比例
%   min_rx_aper   - 最小接收孔径（阵元数），用于限制近场
%   max_rx_aper   - 最大接收孔径（阵元数），用于限制远场
%   win_type      - 窗函数类型: 'hamming', 'hanning', 'chebwin'
%   varargin      - 窗函数的额外参数（如chebwin的旁瓣衰减）
%
% 输出:
%   apod_weights  - 变迹权重表 [nElement, nSample]
%   active_info   - 结构体，包含各采样点的激活信息
min_rx_aper = min_max_aper(1);
max_rx_aper = min_max_aper(2);

% 初始化权重表
apod_weights = zeros(nElement, nSample);

% 中心阵元索引
center_idx = floor((nElement + 1) / 2);

% 存储激活信息
active_aper = zeros(1, nSample);  % 实际使用的孔径
active_start = zeros(1, nSample); % 起始阵元索引
active_end = zeros(1, nSample);   % 结束阵元索引

% 主循环：为每个采样点计算权重
for s = 1:nSample
    % 当前深度对应的有效孔径
    depth = s/fs*c/2;  % 假设采样点索引与深度成正比
    target_aper = round(depth / rx_fnumber/pitch);
    
    % 应用孔径限制
    target_aper = max(target_aper, min_rx_aper);
    target_aper = min(target_aper, max_rx_aper);
    target_aper = min(target_aper, nElement);  % 不能超过总阵元数
    
    % 确保孔径为整数且为奇数（对称）
    aper_actual = floor(target_aper);
    if mod(aper_actual, 2) == 0
        aper_actual = aper_actual - 1;  % 调整为奇数
    end
    aper_actual = max(aper_actual, 1);  % 至少1个阵元
    
    % 计算激活阵元范围
    half_aper = floor(aper_actual / 2);
    start_idx = max(center_idx - half_aper, 1);
    end_idx = min(center_idx + half_aper, nElement);
    
    % 记录激活信息
    active_aper(s) = end_idx - start_idx + 1;
    active_start(s) = start_idx;
    active_end(s) = end_idx;
    
    % 生成窗函数
    switch lower(win_type)
        case 'hamming'
            if isempty(varargin)
                win = hamming(active_aper(s));
            else
                win = hamming(active_aper(s), varargin{1});
            end
            
        case 'hanning'
            if isempty(varargin)
                win = hanning(active_aper(s));
            else
                win = hanning(active_aper(s), varargin{1});
            end
            
        case 'chebwin'
            if ~isempty(varargin)
                % chebwin需要旁瓣衰减参数
                win = chebwin(active_aper(s), varargin{1});
            else
                win = chebwin(active_aper(s), 60);  % 默认60dB衰减
            end
            
        otherwise
            error('不支持的窗函数类型: %s', win_type);
    end
    
    % 应用窗函数
    apod_weights(start_idx:end_idx, s) = win;
end


% % 列归一化：对每个采样点（列）单独归一化
% for s = 1:nSample
%     if active_aper(s) > 0
%         % 只对非零值进行归一化
%         non_zero_mask = apod_weights(:, s) ~= 0;
%         if any(non_zero_mask)
%             non_zero_weights = apod_weights(non_zero_mask, s);
%             weight_sum = sum(non_zero_weights);
%             if weight_sum > 0
%                 % 归一化并保持0值不变
%                 apod_weights(non_zero_mask, s) = non_zero_weights / weight_sum;
%             end
%         end
%     end
% end

% 返回激活信息
if nargout > 1
    active_info.aper = active_aper;
    active_info.start_idx = active_start;
    active_info.end_idx = active_end;
end
end