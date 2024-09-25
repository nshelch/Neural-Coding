% Vibration analysis
load(fullfile(DataPath(), 'ProcessedData', 'Vibrations_data'))
num_electrodes = 256;
cmap = LoadSubjectChannelMap('BCI02');

%% Format spike counts by condition and electrode
% Get unique conditions
uf = unique(data.Freq);
ua = unique(data.Amp);
num_combs = length(uf) * length(ua);
c = lines(num_combs);

% Binning for later plots
bin_size = 0.1;
bin_edges = [-1:bin_size:2];
bin_center = bin_edges(1:end-1) + (bin_size / 2);
analog_time = linspace(-1, 2, 30e3);

% Allocate to cell array of appropriate size
[spike_times, isis] = deal(cell(length(uf), length(ua), num_electrodes));
spike_counts = NaN(size(spike_times));
[vibrometer_data, force_data, accel_data] = deal(cell(length(uf), length(ua)));
for f = 1:length(uf)
    for a = 1:length(ua)
        % Find matching indices
        idx = data.Freq == uf(f) & data.Amp == ua(a);
        if sum(idx) == 0
            continue
        end
        % Spikes
        temp = cat(1, data.Spikes{idx});
        for e = 1:num_electrodes
            spike_times{f,a,e} = temp(:,e);
            temp_isis = cellfun(@(c) pdist(c'), temp(:,e), 'UniformOutput', false);
            isis{f,a,e} = cat(2, temp_isis{:});
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


%% Plot one condition with stimulus
e = 207;
clf;
f = 3;
a = 13;

if ismember(e, cmap.ArrayLocations{1}) || ismember(e, cmap.ArrayLocations{3})
    array_loc = 'Motor';
else
    array_loc = 'Sensory';
end

set(gcf, 'Units', 'Inches', 'Position', [30 1 4 10])
axes('Position', [.1 .75 .8 .2]); hold on
    t = FoldedPSTH(spike_times{f,a,e}, bin_edges, 5, 5);
    AlphaLine(bin_center, t, [.6 .6 .6])
    ylabel('FR')
    title(sprintf('%s: Electrode %d', array_loc, e))

axes('Position', [.1 .4 .8 .3]); hold on
    r = Rasterfy(spike_times{f,a,e}, .4, true);
    plot(r(:,1), r(:,2), 'Color', [.6 .6 .6])
    set(gca, 'YLim', [0 length(spike_times{f,a,e})])

axes('Position', [.1 .05 .8 .3]); hold on
    AlphaLine(analog_time, rescale(force_data{f,a}), rgb(244, 67, 54))
    AlphaLine(analog_time, rescale(accel_data{f,a}) + 1, rgb(103, 58, 183))
    AlphaLine(analog_time, rescale(vibrometer_data{f,a}) + 2, rgb(3, 155, 229))
    set(gca, 'YTick', [])
    text(-0.75, 0, 'Force', 'Color', rgb(244, 67, 54))
    text(-0.75, 1, 'Accelerometer', 'Color', rgb(103, 58, 183))
    text(-0.75, 2, 'Vibrometer', 'Color', rgb(3, 155, 229))
    xlabel('Time (s)')

shg
%% Plot all raster + PSTH
for e = 207%1:num_electrodes
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
                if isempty(spike_times{f,a,e})
                    continue
                end
                t = FoldedPSTH(spike_times{f,a,e}, bin_edges, 5, 5);
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
                if isempty(spike_times{f,a,e})
                    continue
                end
                r = Rasterfy(spike_times{f,a,e}, .4, true);
                plot(r(:,1), r(:,2) + yo, 'Color', c(ci,:))
                yo = yo + length(spike_times{f,a,e});
                ci = ci + 1;
            end
        end
        set(gca, 'YLim', [0 yo])
        xlabel('Time (s)')
    % 
    % print(gcf, fullfile('C:\Users\somlab\Desktop\ExportedImages\NeuralCoding\Vibrations', ...
    %     sprintf('%s_%d', array_loc, e)), '-dpng', '-r300')
end

