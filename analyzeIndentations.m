clear; close all; clc;

addpath(genpath('HelperFunctions'))

%% Set params and load data
load(fullfile(DataPath(), 'ProcessedData', 'Indentations_data'))
numElectrodes = 256;

% Get unique conditions
freq = unique(data.Freq); % mm/s
amp = unique(data.Amp); % mm

%% Format spike counts by condition and electrode

% Allocate to cell array of appropriate size
spikeCounts = cell(length(freq), length(amp), numElectrodes);
[vibrometerData, forceData, accelData] = deal(cell(length(freq), length(amp)));

for f = 1:length(freq)
    for a = 1:length(amp)

        % Find matching indices
        idx = data.Freq == freq(f) & data.Amp == amp(a);

        % Spike data
        tmpSpikes = cat(1, data.Spikes{idx}); % num trials x num electrodes matrix
        for e = 1:numElectrodes
            spikeCounts{f, a, e} = tmpSpikes(:, e);
        end

        % Stimulus data
        tmpStim = cat(1, data.AnalogData(idx));

        % Pad stimulus values
        maxLen = max(cellfun(@length, tmpStim));
        for i = 1:length(tmpStim)
            if length(tmpStim{i}) < maxLen
                padLen = maxLen - length(tmpStim{i});
                tmpStim{i} = cat(2, tmpStim{i}, NaN(3, padLen));
            end
        end

        % Concatenate across trials
        tmpStim = cat(3, tmpStim{:});
        vibrometerData{f, a} = squeeze(tmpStim(1,:,:)); % timepoints x trials
        forceData{f, a} = squeeze(tmpStim(2,:,:));
        accelData{f, a} = squeeze(tmpStim(3,:,:));

    end % amp loop
end % freq loop

%% Compute response to stimulus

baselineInt = [-1 -.5];
onsetInt = [0 0.5];
sustainInt = [.5 .7];
offsetInt = [1 1.2];

% stimModValue = zeros(length(freq), length(amp), numElectrodes);
% adaptationIndex = zeros(length(freq), length(amp), numElectrodes);

% for f = 1:length(freq)
%     for a = 1:length(amp)
%         for e = 1:numElectrodes
%             tmpSpikes = spikeCounts{f, a, e};
%
%             baselineSC = cellfun(@(z) sum(z(:) >= baselineInt(1) & z(:) <= baselineInt(2)), tmpSpikes);
%             onsetSC = cellfun(@(z) sum(z(:) >= onsetInt(1) & z(:) <= onsetInt(2)), tmpSpikes);
%             sustainSC = cellfun(@(z) sum(z(:) >= sustainInt(1) & z(:) <= sustainInt(2)), tmpSpikes);
%             offsetSC = cellfun(@(z) sum(z(:) >= offsetInt(1) & z(:) <= offsetInt(2)), tmpSpikes);
%
%             baselineFR = baselineSC ./ diff(baselineInt);
%             onsetFR = onsetSC ./ diff(onsetInt);
%             sustainFR = sustainSC ./ diff(sustainInt);
%             offsetFR = offsetSC ./ diff(offsetInt);
%
% %             stimModValue(f, a, e) = (mean(onsetFR) - mean(baselineFR)) ./ std(baselineFR);
%
%             stimModValue(f, a, e) = (median(onsetFR) - median(baselineFR)) ./ max(abs(baselineFR - median(baselineFR)));
%
%             adaptationIndex(f, a, e) = computeAdaptationIndex(median(offsetFR - median(baselineFR)), median(sustainFR - median(baselineFR)), 'norm');
%
%         end % electrode loop
%     end % amp loop
% end % freq loop

stimModValue = zeros(1, numElectrodes);
adaptationIndex = zeros(1, numElectrodes);

for e = 1:numElectrodes
    tmpSpikes = [spikeCounts{:, :, e}];

    baselineSC = cellfun(@(z) sum(z(:) >= baselineInt(1) & z(:) <= baselineInt(2)), tmpSpikes);
    onsetSC = cellfun(@(z) sum(z(:) >= onsetInt(1) & z(:) <= onsetInt(2)), tmpSpikes);
    sustainSC = cellfun(@(z) sum(z(:) >= sustainInt(1) & z(:) <= sustainInt(2)), tmpSpikes);
    offsetSC = cellfun(@(z) sum(z(:) >= offsetInt(1) & z(:) <= offsetInt(2)), tmpSpikes);

    baselineFR = baselineSC ./ diff(baselineInt);
    onsetFR = onsetSC ./ diff(onsetInt);
    sustainFR = sustainSC ./ diff(sustainInt);
    offsetFR = offsetSC ./ diff(offsetInt);

    baselineMu = mean(baselineFR(:));
    baselineSigma = std(baselineFR(:));
    stimModValue(e) = (mean(onsetFR(:)) - baselineMu) ./ baselineSigma;

    %         baselineMu = median(baselineFR(:));
    %     baselineSigma = max(abs(baselineFR(:) - baselineMu));
    %     stimModValue(e) = (median(onsetFR) - baselineMu) ./ baselineSigma;

    adaptationIndex(e) = computeAdaptationIndex(mean(offsetFR(:)) - baselineMu, mean(sustainFR(:)) - baselineMu, 'tan');

end % electrode loop

%% Saving data

save('./Data/HumanIndentSCData.mat', )


%% Summary plots

modThresh = 1.5;

sensoryElecIdx = [65:96, 193:224];

histogram(adaptationIndex(sensoryElecIdx(stimModValue(sensoryElecIdx) >= modThresh )), 'BinEdges', linspace(0, 1, 7))

%% Plotting

chanMap = LoadSubjectChannelMap('BCI02');

% Binning for later plots
dt = 0.1;
binEdges = -1:dt:2;
binCenter = binEdges(1:end-1) + (dt / 2);

analogTime = linspace(-1, 2, 30e3);

% Plot one condition with stimulus
for e = 71
    clf;
    f = 5;
    a = 5;

    if ismember(e, chanMap.ArrayLocations{1}) || ismember(e, chanMap.ArrayLocations{3})
        array_loc = 'Motor';
        plt_color = rgb(0, 91, 127);
    else
        array_loc = 'Sensory';
        plt_color = rgb(152, 3, 126);
    end

    set(gcf, 'Units', 'Inches', 'Position', [10 1 4 9])
    axes('Position', [.15 .75 .8 .2]); hold on
    t = FoldedPSTH(spikeCounts{f, a, e}, binEdges, 5, 5);
    AlphaLine(binCenter, t, plt_color)
    ylabel('FR')
    title(sprintf('%s: Electrode %d', array_loc, e))
    yLim = get(gca, 'YLim');

    if exist('stimModValue', 'var')
        text(-1, yLim(2), {sprintf('MV: %.2f', stimModValue(f, a, e)); sprintf('AI: %.2f', adaptationIndex(f, a, e))}, ...
            'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    end

    axes('Position', [.15 .4 .8 .3]); hold on
    r = Rasterfy(spikeCounts{f, a, e}, .4, true);
    plot(r(:,1), r(:,2), 'Color', plt_color)
    set(gca, 'YLim', [0 length(spikeCounts{f,a,e})], 'YTick', [])
    xlabel('Time (seconds)')

    axes('Position', [.15 .05 .8 .3]); hold on
    AlphaLine(analogTime, rescale(forceData{f, a}), rgb(244, 67, 54))
    AlphaLine(analogTime, rescale(accelData{f, a}) + 1, rgb(103, 58, 183))
    AlphaLine(analogTime, rescale(vibrometerData{f, a}) + 2, rgb(3, 155, 229))
    set(gca, 'YTick', [])
    text(-0.75, 0, 'Force', 'Color', rgb(244, 67, 54))
    text(-0.75, 1, 'Accelerometer', 'Color', rgb(103, 58, 183))
    text(-0.75, 2, 'Vibrometer', 'Color', rgb(3, 155, 229))
    shg

    %     print(gcf, fullfile('.\Figures\Indentations', ...
    %         sprintf('%s_%d', array_loc, e)), '-dpng', '-r300')
end
%% Plot all raster + PSTH

c = lines(25);

for e = 216 %1:numElectrodes
    if ismember(e, chanMap.ArrayLocations{1}) || ismember(e, chanMap.ArrayLocations{3})
        array_loc = 'Motor';
    else
        array_loc = 'Sensory';
    end

    clf;
    set(gcf, 'Units', 'Inches', 'Position', [10 1 4 9])
    axes('Position', [.1 .75 .8 .2]); hold on
    ci = 1;
    for f = 1:length(freq)
        for a = 1:length(amp)
            t = FoldedPSTH(spikeCounts{f,a,e}, binEdges, 5, 5);
            AlphaLine(binCenter, t, c(ci,:))
            ci = ci + 1;
        end
    end
    ylabel('FR')
    title(sprintf('%s: Electrode %d', array_loc, e))

        yLim = get(gca, 'YLim');
    if exist('stimModValue', 'var')
        text(-.95, yLim(2), {sprintf('MV: %.2f', stimModValue(e)); sprintf('AI: %.2f', adaptationIndex(e))}, ...
            'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    end

    axes('Position', [.1 .1 .8 .6]); hold on
    yo = 0; ci = 1;
    for f = 1:length(freq)
        for a = 1:length(amp)
            r = Rasterfy(spikeCounts{f,a,e}, .4, true);
            plot(r(:,1), r(:,2) + yo, 'Color', c(ci,:))
            yo = yo + length(spikeCounts{f,a,e});
            ci = ci + 1;
        end
    end
    set(gca, 'YLim', [0 yo])
    xlabel('Time (s)')

%     print(gcf, fullfile('.\Figures\Indentations', ...
%         sprintf('%s_%d', array_loc, e)), '-dpng', '-r300')

end

