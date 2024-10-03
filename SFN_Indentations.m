%% Human Indentation Data
addpath(genpath('HelperFunctions'))

% Set params and load data
load(fullfile(DataPath(), 'ProcessedData', 'Indentations_data'))
numElectrodes = 256;

% Get unique conditions
freq = unique(data.Freq); % mm/s
amp = unique(data.Amp); % mm

% Format spike counts by condition and electrode

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

% Compute response to stimulus

baselineInt = [-1 -.5];
onsetInt = [0 0.5];
sustainInt = [.5 .7];
offsetInt = [1 1.2];

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

    adaptationIndex(e) = computeAdaptationIndex(mean(offsetFR(:)) - baselineMu, mean(sustainFR(:)) - baselineMu, 'tan');

end % electrode loop

%% Example Rasters - Human

% Binning for later plots
dt = 0.1;
binEdges = -1:dt:2;
binCenter = binEdges(1:end-1) + (dt / 2);

analogTime = linspace(-1, 2, 30e3);

humanColor = rgb(67, 160, 71);

%% Plot one condition with stimulus
% close all; clc;
coi = [71, 207, 214, 216, 218]; 

figure('Units', 'Inches', 'Position', [10 2 4 6])
for e = coi
    clf;
    f = 5;
    a = 5;

    axes('Position', [.15 .75 .8 .2]); hold on
    t = FoldedPSTH(spikeCounts{f, a, e}, binEdges, 5, 5);
    AlphaLine(binCenter, t, humanColor)
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', [])

    axes('Position', [.15 .1 .8 .6]); hold on
    r = Rasterfy(spikeCounts{f, a, e}, .4, true);
    plot(r(:,1), r(:,2), 'Color', humanColor)
    set(gca, 'YLim', [0 length(spikeCounts{f,a,e})], 'YTick', [], 'YColor', 'none')
    xlabel('Time (seconds)')
    print(gcf, fullfile('.\Figures\SFN', sprintf('HumanIndent_SE%d', e)), '-dpng', '-r900')

end


%% Plot all raster + PSTH
coi = [195]; % inhibitory response: 93, 212, 214, 216, 218

c = cmap_gradient([rgb(76, 175, 80); rgb(27, 94, 32)], 5);

close all;
figure('Units', 'Inches', 'Position', [10 2 4 6])
for e = 195
    clf;    
    axes('Position', [.15 .75 .8 .2]); hold on
    ci = 1;
    for f = 5 %1:length(freq)
        for a = 1:length(amp)
            t = FoldedPSTH(spikeCounts{f,a,e}, binEdges, 5, 5);
            AlphaLine(binCenter, t, c(ci,:))
            ci = ci + 1;
        end
    end
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', [])

    axes('Position', [.15 .1 .8 .6]); hold on
    yo = 0; ci = 1;
    for f = 5 %1:length(freq)
        for a = 1:length(amp)
            r = Rasterfy(spikeCounts{f,a,e}, .4, true);
            plot(r(:,1), r(:,2) + yo, 'Color', c(ci,:))
            yo = yo + length(spikeCounts{f,a,e});
            ci = ci + 1;
        end
    end
    set(gca, 'YLim', [0 yo], 'YTick', [], 'YColor', 'none')
    xlabel('Time (s)')

%     print(gcf, fullfile('.\Figures\SFN', sprintf('HumanIndent_SE%d', e)), '-dpng', '-r900')

end

%% Plot averaged responses

 % inhibitory response: 214, 216, 218
humanColor = rgb(67, 160, 71);

close all;
figure('Units', 'Inches', 'Position', [8 4 3 2])
for e = 212
    clf;    
    axes('Position', [.175 .25 .8 .725]); hold on
    for f = 5 
            t = FoldedPSTH(vertcat(spikeCounts{f, 1:5, e}), binEdges, 5, 5);
            AlphaLine(binCenter, t, humanColor)
    end
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', -1:1:2)
    xlabel('Time (s)')
    
%     % Stimulus lines
%     plot([0 0], [25 40], '--', 'Color', [200, 200, 200] ./ 255)
%     plot([1 1], [25 40], '--', 'Color', [200, 200, 200] ./ 255)

    print(gcf, fullfile('.\Figures\SFN', sprintf('HumanIndent_SE%d', e)), '-dpng', '-r900')
end

%%
close all;
figure('Units', 'Inches', 'Position', [8 4 3 4])
for e = 212
    clf;    
    axes('Position', [.175 .125 .8 .525]); hold on
    for f = 5 
            t = FoldedPSTH(vertcat(spikeCounts{f, 1:5, e}), binEdges, 5, 5);
            AlphaLine(binCenter, t, humanColor)
    end
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', -1:1:2)
    xlabel('Time (s)')
    
%     % Stimulus lines
    axes('Position', [.175 .75 .8 .2]); hold on
    AlphaLine(analogTime, rescale(forceData{5, 5}), rgb(0, 0, 0), ...
        'ErrorAlpha', 0, 'EdgeAlpha', 0, 'LineWidth', 1.5)
    set(gca, 'YTick', [], 'YColor', 'none', 'XTickLabel', '', 'XColor', 'none')
    ylabel('Force', 'Color', 'k')
    text(-.95, 1, 'Stimulus', ...
        'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

%     print(gcf, fullfile('.\Figures\SFN', sprintf('HumanIndent_SE%d_V3', e)), '-dpng', '-r900')
end



%% Plot human motor data

motorColor = rgb(0, 137, 123);

close all;
figure('Units', 'Inches', 'Position', [8 4 5 4])
for e = 154
    clf;    
    axes('Position', [.15 .15 .8 .8]); hold on
    for f = 5 %1:length(freq)
            t = FoldedPSTH(vertcat(spikeCounts{f, 1:5, e}), binEdges, 5, 5);
            AlphaLine(binCenter, t, motorColor)
    end
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', -1:1:2)
    xlabel('Time (s)')
    
%     yLim = get(gca, 'YLim');
%     patch([0 1 1 0], [yLim(1) yLim(1) yLim(2) yLim(2)], [200, 200, 200] ./ 255, 'EdgeColor', 'none')



%     print(gcf, fullfile('.\Figures\SFN', sprintf('HumanIndent_ME%d', e)), '-dpng', '-r900')

end




%% Load macaque data

load('S:\UserFolders\CharlesGreenspon\SpeciesVibrations\RawMacaqueCorticalIndentationData.mat')

a1Idx = find(strcmp(L4_str.specs.anat_loc, 'a1'));
macaqueAI = L4_str.data.transienceIndex(a1Idx);

nhpColor = rgb(251, 140, 0);

%%
% Example Rasters - Macaque
figure('Units', 'Inches', 'Position', [10 2 4 6])

coi = [14, 72, 117];
for e = coi
    clf;
    axes('Position', [.15 .75 .8 .2]); hold on
    t = FoldedPSTH(cellfun(@transpose, L4_str.data.rasters{e}, 'UniformOutput', false), binEdges, 5, 5);
    AlphaLine(binCenter, t, nhpColor)
    ylabel('Firing Rate (Hz)')
    set(gca, 'XLim', [-1 2], 'XTick', [])

    axes('Position', [.15 .1 .8 .6]); hold on
    r = Rasterfy(cellfun(@transpose, L4_str.data.rasters{e}, 'UniformOutput', false), .4, true);
    plot(r(:,1), r(:,2), 'Color', nhpColor)
    set(gca, 'XLim', [-1 2], 'YLim', [0 length(r)], 'YTick', [], 'YColor', 'none')
    xlabel('Time (seconds)')
    
    print(gcf, fullfile('.\Figures\SFN', sprintf('MacaqueIndent_%d', e)), '-dpng', '-r900')

end

figure('Units', 'Inches', 'Position', [10 2 4 6])
r = Rasterfy(indentation_data.aligned_spikes(1, :), .4, true);
plot(r(:,1), r(:,2), 'Color', nhpColor)


%% Adaptibility Index Figure

nhpColor = rgb(251, 140, 0);
humanColor = rgb(67, 160, 71);

sensoryChIdx = [65:96, 193:224];

figure('Units', 'Inches', 'Position', [10 2 4 5])

axes('Position', [.15 .55 .8 .375])
histogram(adaptationIndex(sensoryChIdx), 'BinEdges', linspace(0, 1, 7), 'Normalization', 'probability', ...
    'FaceColor', humanColor)
box off
set(gca, 'YTick', 0:.2:.6, 'XTick', [])
text(0, .6, ColorText({'Human - Area 1'}, humanColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')


axes('Position', [.15 .1 .8 .375])
histogram(macaqueAI, 'BinEdges', linspace(0, 1, 7), 'Normalization', 'probability', ...
    'FaceColor', nhpColor)
box off
set(gca, 'YTick', 0:.1:.3, 'XTick', [0 1])
text(0, .3, ColorText({'NHP - Area 1'}, nhpColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
xlabel('Adaptation Index')
ylabel('Proportion', 'Position', [-.2 .325 -1])
   
print(gcf, fullfile('.\Figures\SFN', 'adaptationIndex'), '-dpng', '-r900')



histogram(indentation_data.adaptability_index, 'BinEdges', linspace(0, 1, 7), 'Normalization', 'probability', ...
    'FaceColor', nhpColor)

