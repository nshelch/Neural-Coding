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

baselineInt = [-1, -.6]; % [-1 -.5];
baseline2Int = [-.6 -.2];
onsetInt = [0 0.4]; % [0 0.5];
sustainInt = [.5 .7];
offsetInt = [1 1.2];

% Permutation Params
numReps = 1000;
numSamples = 100;
alpha = 0.1;

nullDist = zeros(numElectrodes, numReps);
pMod = zeros(1, numElectrodes);
isMod = zeros(1, numElectrodes);
stimModValue = zeros(1, numElectrodes);
adaptationIndex = zeros(1, numElectrodes);

for e = 1:numElectrodes
    tmpSpikes = [spikeCounts{:, :, e}];

    baselineSC = cellfun(@(z) sum(z(:) >= baselineInt(1) & z(:) <= baselineInt(2)), tmpSpikes);
    baseline2SC = cellfun(@(z) sum(z(:) >= baseline2Int(1) & z(:) <= baseline2Int(2)), tmpSpikes);

    onsetSC = cellfun(@(z) sum(z(:) >= onsetInt(1) & z(:) <= onsetInt(2)), tmpSpikes);
    sustainSC = cellfun(@(z) sum(z(:) >= sustainInt(1) & z(:) <= sustainInt(2)), tmpSpikes);
    offsetSC = cellfun(@(z) sum(z(:) >= offsetInt(1) & z(:) <= offsetInt(2)), tmpSpikes);

    baselineFR = baselineSC ./ diff(baselineInt);
    baseline2FR = baseline2SC ./ diff(baseline2Int);
    onsetFR = onsetSC ./ diff(onsetInt);
    sustainFR = sustainSC ./ diff(sustainInt);
    offsetFR = offsetSC ./ diff(offsetInt);

    baselineMu(e) = mean(baselineFR(:));
    baselineSigma(e) = std(baselineFR(:));

    % Permutation Test
    baseline2FR = baseline2FR(:);
    for r = 1:numReps
        tmp = datasample(baseline2FR, numSamples, 'Replace', true);
        nullDist(e, r) = (mean(tmp) - baselineMu(e)) ./ baselineSigma(e);
    end % permutation loop

    % Significance test?
    stimModValue(e) = (mean(onsetFR(:)) - baselineMu(e)) ./ baselineSigma(e);

    leftTail = 1 - sum(stimModValue(e) < nullDist(e, :)) / numReps;
    rightTail = 1 - (sum(stimModValue(e) > nullDist(e, :)) / numReps);
    pMod(e) = min([leftTail, rightTail]);
    isMod(e) = pMod(e) <= alpha / 2;

    %         baselineMu = median(baselineFR(:));
    %     baselineSigma = max(abs(baselineFR(:) - baselineMu));
    %     stimModValue(e) = (median(onsetFR) - baselineMu) ./ baselineSigma;

    % Adaptation Index
    adaptationIndex(e) = computeAdaptationIndex(mean(offsetFR(:)) - baselineMu(e), mean(sustainFR(:)) - baselineMu(e), 'tan');

end % electrode loop

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


%% CCR Analysis

% Get the significantly modulated channels

sigElecIdx = find(baselineMu >= 5);
numSigElecs = length(sigElecIdx);

% Find the sensory and motor channel indices
chanMap = LoadSubjectChannelMap('BCI02');

chFilter = cat(1, chanMap.ArrayLocations{chanMap.IsSensory});
sensoryChannelInd = sort(chFilter(~isnan(chFilter)));

chFilter = cat(1, chanMap.ArrayLocations{chanMap.IsMotor});
motorChannelInd = sort(chFilter(~isnan(chFilter)));

% Order the significantly modulated electrodes into sensory and motor
% groups
elecOrder = [sigElecIdx(ismember(sigElecIdx, sensoryChannelInd)), sigElecIdx(ismember(sigElecIdx, motorChannelInd))];

% Binning params
dt = 0.01;
binEdges = [-.2:dt:1.5];
gaussWidth = 5; % for smoothing
numTrials = 1500; % Need to compute from data

maxCorr = NaN(numSigElecs, numSigElecs);
maxLag = NaN(numSigElecs, numSigElecs);

for n1 = 1:numSigElecs
    tmpSpikes = vertcat(spikeCounts{1:5, 1:5, elecOrder(n1)});
    tmpSpikes = horzcat(tmpSpikes{:});
    n1PSTH = (histcounts(tmpSpikes, binEdges) ./ dt) ./ numTrials ;
    n1PSTH = (n1PSTH - baselineMu(elecOrder(n1))) ./ baselineSigma(elecOrder(n1));

    for n2 = (n1 + 1):numSigElecs
        tmpSpikes = vertcat(spikeCounts{1:5, 1:5, elecOrder(n2)});
        tmpSpikes = horzcat(tmpSpikes{:});
        n2PSTH = (histcounts(tmpSpikes, binEdges) ./ dt) ./ numTrials ;
        n2PSTH = (n2PSTH - baselineMu(elecOrder(n2))) ./ baselineSigma(elecOrder(n2));

        % Compute the CCR and get the max corr and lag which corresponds to
        % it
        [corrVal, lags] = xcorr(n1PSTH, n2PSTH, 10, 'normalized');
        [~, lagIdx] = max(abs(corrVal));
        maxCorr(n1, n2) = corrVal(lagIdx); maxCorr(n2, n1) = corrVal(lagIdx);
        maxLag(n1, n2) = lags(lagIdx);     maxLag(n2, n1) = lags(lagIdx);

        % Check individual channel results
%         subplot(3,1,1)
%         plot(binEdges(1:end - 1) + (dt / 2), n1PSTH)
%         subplot(3,1,2)
%         plot(binEdges(1:end - 1) + (dt / 2), n2PSTH)
%         subplot(3,1,3)
%         stem(lags, corrVal)

    end % channel 2 loops

end % channel 1 loop

%%

cmap = cmap_gradient([rgb(25, 118, 210); rgb(250, 250, 250); rgb(194, 24, 91)]);

maxCorr2 = maxCorr;
maxCorr2(isnan(maxCorr)) = 1;

figure;
h = imagesc(maxCorr2);
set(h, 'AlphaData', ~isnan(maxCorr2))
axis square
cbar = colorbar;
ylabel(cbar, 'Correlation')
cbar.Ticks = -1:2:1;
cbar.TickLabels = {'', '', ''};
cbar.Limits = [-1 1];


hold on

% Anterior sensory array line
chFilter = cat(1, chanMap.ArrayLocations{2});
chFilter = sort(chFilter(~isnan(chFilter)));
antSensoryChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');

plot([0.5, length(elecOrder) + 0.5], [antSensoryChIdx + 0.5, antSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
plot([antSensoryChIdx + 0.5, antSensoryChIdx + 0.5], [0.5 length(elecOrder) + 0.5], 'k', 'LineWidth', 1) % right

% Posterior sensory array line
chFilter = cat(1, chanMap.ArrayLocations{4});
chFilter = sort(chFilter(~isnan(chFilter)));
postSensoryChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');

plot([0.5, length(elecOrder) + 0.5], [postSensoryChIdx + 0.5, postSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
plot([postSensoryChIdx + 0.5, postSensoryChIdx + 0.5], [0.5 length(elecOrder) + 0.5], 'k', 'LineWidth', 1) % right

% Anterior motor array line
chFilter = cat(1, chanMap.ArrayLocations{1});
chFilter = sort(chFilter(~isnan(chFilter)));
antMotorChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');

plot([0.5, length(elecOrder) + 0.5], [antMotorChIdx + 0.5, antMotorChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
plot([antMotorChIdx + 0.5, antMotorChIdx + 0.5], [0.5 length(elecOrder) + 0.5], 'k', 'LineWidth', 1) % right

% Posterior motor array line
chFilter = cat(1, chanMap.ArrayLocations{3});
chFilter = sort(chFilter(~isnan(chFilter)));
postMotorChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');

plot([0.5, length(elecOrder) + 0.5], [postMotorChIdx + 0.5, postMotorChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
plot([postMotorChIdx + 0.5, postMotorChIdx + 0.5], [0.5 length(elecOrder) + 0.5], 'k', 'LineWidth', 1) % right


% % Anterior sensory array line
% chFilter = cat(1, chanMap.ArrayLocations{2});
% chFilter = sort(chFilter(~isnan(chFilter)));
% antSensoryChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');
% 
% plot([0.5, antSensoryChIdx + 0.5], [antSensoryChIdx + 0.5, antSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
% plot([antSensoryChIdx + 0.5, antSensoryChIdx + 0.5], [antSensoryChIdx + 0.5, 0], 'k', 'LineWidth', 1) % right
% 
% % Posterior sensory array line
% chFilter = cat(1, chanMap.ArrayLocations{4});
% chFilter = sort(chFilter(~isnan(chFilter)));
% postSensoryChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');
% 
% plot([antSensoryChIdx + 0.5, postSensoryChIdx + 0.5], [antSensoryChIdx + 0.5, antSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % top
% plot([antSensoryChIdx + 0.5, antSensoryChIdx + 0.5], [postSensoryChIdx + 0.5, antSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % left
% plot([antSensoryChIdx + 0.5, postSensoryChIdx + 0.5], [postSensoryChIdx + 0.5, postSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
% plot([postSensoryChIdx + 0.5, postSensoryChIdx + 0.5], [postSensoryChIdx + 0.5, antSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % right
% 
% % Anterior motor array line
% chFilter = cat(1, chanMap.ArrayLocations{1});
% chFilter = sort(chFilter(~isnan(chFilter)));
% antMotorChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');
% 
% plot([postSensoryChIdx + 0.5, antMotorChIdx + 0.5], [postSensoryChIdx + 0.5, postSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % top
% plot([postSensoryChIdx + 0.5, postSensoryChIdx + 0.5], [antMotorChIdx + 0.5, postSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % left
% plot([postSensoryChIdx + 0.5, antMotorChIdx + 0.5], [antMotorChIdx + 0.5, antMotorChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
% plot([antMotorChIdx + 0.5, antMotorChIdx + 0.5], [antMotorChIdx + 0.5, postSensoryChIdx + 0.5], 'k', 'LineWidth', 1) % right
% 
% % Posterior motor array line
% chFilter = cat(1, chanMap.ArrayLocations{3});
% chFilter = sort(chFilter(~isnan(chFilter)));
% postMotorChIdx = find(ismember(elecOrder,  chFilter), 1, 'last');
% 
% plot([antMotorChIdx + 0.5, postMotorChIdx + 0.5], [antMotorChIdx + 0.5, antMotorChIdx + 0.5], 'k', 'LineWidth', 1) % top
% plot([antMotorChIdx + 0.5, antMotorChIdx + 0.5], [postMotorChIdx + 0.5, antMotorChIdx + 0.5], 'k', 'LineWidth', 1) % left
% plot([antMotorChIdx + 0.5, postMotorChIdx + 0.5], [postMotorChIdx + 0.5, postMotorChIdx + 0.5], 'k', 'LineWidth', 1) % bottom
% plot([postMotorChIdx + 0.5, postMotorChIdx + 0.5], [postMotorChIdx + 0.5, antMotorChIdx + 0.5], 'k', 'LineWidth', 1) % right

set(gca, 'YTick', [7, 20, 32, 65], 'YTickLabel', {'Ant. SC', 'Post. SC', 'Ant. MC', 'Post. MC'})
set(gca, 'XTick', [7, 20, 32, 65], 'XTickLabel', {'Ant. SC', 'Post. SC', 'Ant. MC', 'Post. MC'})
set(gca, 'XTick', [], 'XTickLabel', {''})
title('Max Correlation')

% nexttile
% h = imagesc(maxLag);
% set(h, 'AlphaData', ~isnan(maxLag))
% axis square
% cbar = colorbar;
% ylabel(cbar, 'Lag Time')
% set(gca, 'YTick', [1, 5, 30, 58], 'YTickLabel', {'Ant. SC', 'Post. SC', 'Ant. MC', 'Post. MC'})
% title('Lag Time')


%% Saving data

% save('./Data/HumanIndentSCData.mat', )


%% Plotting

chanMap = LoadSubjectChannelMap('BCI02');

% Binning for later plots
dt = 0.1;
binEdges = -1:dt:2;
binCenter = binEdges(1:end-1) + (dt / 2);

analogTime = linspace(-1, 2, 30e3);

%% Plot Average ACROSS Condition Response

figure('Units', 'Inches', 'Position', [8 4 4 3])
for e = 1:256
    clf;

    if ismember(e, chanMap.ArrayLocations{1}) || ismember(e, chanMap.ArrayLocations{3})
        array_loc = 'Motor';
        plt_color = rgb(0, 91, 127);
    else
        array_loc = 'Sensory';
        plt_color = rgb(152, 3, 126);
    end

    axes('Position', [.15 .15 .8 .725]); hold on

    t = FoldedPSTH(vertcat(spikeCounts{1:5, 1:5, e}), binEdges, 5, 5);
    AlphaLine(binCenter, t, plt_color)
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', -1:1:2)
    xlabel('Time (s)')

    title(sprintf('%s: Electrode %d', array_loc, e))
    yLim = get(gca, 'YLim');

    if exist('stimModValue', 'var')
        if isMod(e)
            text(-.9, yLim(2), {sprintf('MV: %.2f', stimModValue(e)); sprintf('p: %.2f', pMod(e)); sprintf('AI: %.2f', adaptationIndex(e))}, ...
                'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
        else
            text(-.9, yLim(2), {sprintf('MV: %.2f', stimModValue(e)); sprintf('p: %.2f', pMod(e)); sprintf('AI: %.2f', adaptationIndex(e))}, ...
                'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
        end
    end

    print(gcf, fullfile('.\Figures\Indentations', ...
        sprintf('%s_%d', array_loc, e)), '-dpng', '-r300')

end

%%



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

