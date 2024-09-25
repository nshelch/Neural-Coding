% Indentation analysis
addpath(genpath('HelperFunctions'))
load(fullfile(DataPath(), 'ProcessedData', 'Indentations_data'))
num_electrodes = 256;
cmap = LoadSubjectChannelMap('BCI02');

%% Format spike counts by condition and electrode
% Get unique conditions
uf = unique(data.Freq); % mm/s
ua = unique(data.Amp); % mm
num_combs = length(uf) * length(ua);
c = lines(num_combs);

% Binning for later plots
bin_size = 0.1;
bin_edges = [-1:bin_size:2];
bin_center = bin_edges(1:end-1) + (bin_size / 2);
analog_time = linspace(-1, 2, 30e3);

% Allocate to cell array of appropriate size
spike_counts = cell(length(uf), length(ua), num_electrodes);
[vibrometer_data, force_data, accel_data] = deal(cell(length(uf), length(ua)));
for f = 1:length(uf)
    for a = 1:length(ua)
        % Find matching indices
        idx = data.Freq == uf(f) & data.Amp == ua(a);
        % Spikes
        temp = cat(1, data.Spikes{idx});
        for e = 1:num_electrodes
            spike_counts{f,a,e} = temp(:,e);
        end
        % Stim data
        temp = cat(1, data.AnalogData(idx));
        % pad values
        max_len = max(cellfun(@length, temp));
        for i = 1:length(temp)
            if length(temp{i}) < max_len
                temp{i} = cat(2, temp{i}, NaN(3, 1));
            end
        end
        % Concatenate
        temp = cat(3, temp{:});
        vibrometer_data{f,a} = squeeze(temp(1,:,:));
        force_data{f,a} = squeeze(temp(2,:,:));
        accel_data{f,a} = squeeze(temp(3,:,:));
    end
end

return

%% Plot one condition with stimulus
for e = 71%[71, 214, 216, 232, 240, 245, 233]
    clf;
    f = 5;
    a = 5;
    
    if ismember(e, cmap.ArrayLocations{1}) || ismember(e, cmap.ArrayLocations{3})
        array_loc = 'Motor';
        plt_color = rgb(0, 91, 127);
    else
        array_loc = 'Sensory';
        plt_color = rgb(152, 3, 126);    
    end
    
    set(gcf, 'Units', 'Inches', 'Position', [30 1 4 10])
    axes('Position', [.15 .75 .8 .2]); hold on
        t = FoldedPSTH(spike_counts{f,a,e}, bin_edges, 5, 5);
        AlphaLine(bin_center, t, plt_color)
        ylabel('FR')
        title(sprintf('%s: Electrode %d', array_loc, e))
    
    axes('Position', [.15 .4 .8 .3]); hold on
        r = Rasterfy(spike_counts{f,a,e}, .4, true);
        plot(r(:,1), r(:,2), 'Color', plt_color)
        set(gca, 'YLim', [0 length(spike_counts{f,a,e})], 'YTick', [])
        xlabel('Time (seconds)')
    
    axes('Position', [.15 .05 .8 .3]); hold on
        AlphaLine(analog_time, rescale(force_data{f,a}), rgb(244, 67, 54))
        AlphaLine(analog_time, rescale(accel_data{f,a}) + 1, rgb(103, 58, 183))
        AlphaLine(analog_time, rescale(vibrometer_data{f,a}) + 2, rgb(3, 155, 229))
        set(gca, 'YTick', [])
        text(-0.75, 0, 'Force', 'Color', rgb(244, 67, 54))
        text(-0.75, 1, 'Accelerometer', 'Color', rgb(103, 58, 183))
        text(-0.75, 2, 'Vibrometer', 'Color', rgb(3, 155, 229))
    % 
    % print(gcf, fullfile('C:\Users\somlab\Desktop\ExportedImages\NeuralCoding\Temp', ...
    %         sprintf('%s_%d', array_loc, e)), '-dpng', '-r300')
    
    shg
end
%% Plot all raster + PSTH
for e = 71%1:num_electrodes
    if ismember(e, cmap.ArrayLocations{1}) || ismember(e, cmap.ArrayLocations{3})
        array_loc = 'Motor';
    else
        array_loc = 'Sensory';
    end
    
    clf;
    set(gcf, 'Units', 'Inches', 'Position', [30 1 4 10])
    axes('Position', [.1 .75 .8 .2]); hold on
        ci = 1;
        for f = 1:length(uf)
            for a = 1:length(ua)
                t = FoldedPSTH(spike_counts{f,a,e}, bin_edges, 5, 5);
                AlphaLine(bin_center, t, c(ci,:))
                ci = ci + 1;
            end
        end
        ylabel('FR')
        title(sprintf('%s: Electrode %d', array_loc, e))
    
    axes('Position', [.1 .1 .8 .6]); hold on
        yo = 0; ci = 1;
        for f = 1:length(uf)
            for a = 1:length(ua)
                r = Rasterfy(spike_counts{f,a,e}, .4, true);
                plot(r(:,1), r(:,2) + yo, 'Color', c(ci,:))
                yo = yo + length(spike_counts{f,a,e});
                ci = ci + 1;
            end
        end
        set(gca, 'YLim', [0 yo])
        xlabel('Time (s)')
    
    % print(gcf, fullfile('C:\Users\somlab\Desktop\ExportedImages\NeuralCoding\Indentations', ...
        % sprintf('%s_%d', array_loc, e)), '-dpng', '-r300')
end


%% Plot all conditions
for e = 71%[71, 214, 216, 232, 240, 245, 233]
    clf;
    
    if ismember(e, cmap.ArrayLocations{1}) || ismember(e, cmap.ArrayLocations{3})
        array_loc = 'Motor';
        plt_color = rgb(0, 91, 127);
    else
        array_loc = 'Sensory';
        plt_color = rgb(152, 3, 126);    
    end
    
    set(gcf, 'Units', 'Inches', 'Position', [30 1 4 5])
    axes('Position', [.15 .55 .8 .35]); hold on
        t = FoldedPSTH(cat(2, spike_counts{:,:,e}), bin_edges, 5, 5);
        AlphaLine(bin_center, t, plt_color)
        ylabel('FR')
        title(sprintf('%s: Electrode %d', array_loc, e))
        set(gca, 'XTickLabels', {})
    
    axes('Position', [.15 .1 .8 .4]); hold on
        r = Rasterfy(cat(2, spike_counts{:,:,e}), .4, true);
        plot(r(:,1), r(:,2), 'Color', plt_color)
        set(gca, 'YLim', [0 length(spike_counts{f,a,e})], 'YTick', [])
        xlabel('Time (seconds)')
    

    
    % print(gcf, fullfile('C:\Users\somlab\Desktop\ExportedImages\NeuralCoding\Temp', ...
    %         sprintf('%s_%d', array_loc, e)), '-dpng', '-r300')
    
    shg
end