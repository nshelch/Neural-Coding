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

%% Plot vibration and indentation stimuli

% Vibrations
amplitude = 0.1;
period = 1.0;
t = linspace(0, 10, 5000);
y = amplitude * sin(2 * pi * t / period);
plot(t, y, 'k-', 'LineWidth', 2);
axis off
set(gcf, 'color', 'none'); set(gca, 'color', 'none');
    set(gcf,'InvertHardcopy','off');

print(gcf, fullfile('.\Figures\SFN', 'Correlation'), '-dpng', '-r900')
exportgraphics(gca, fullfile('.\Figures\SFN', 'VibrationStimuli2.png'),'BackgroundColor','none')


%% Example Rasters - Human

% Binning for later plots
dt = 0.1;
binEdges = -1:dt:2;
binCenter = binEdges(1:end-1) + (dt / 2);

analogTime = linspace(-1, 2, 30e3);

humanColor = rgb(57, 73, 171);

%% Plot one condition with stimulus
% close all; clc;
coi = [71, 207, 214, 216, 218]; 

SetFont('Arial', 14)
% figure('Units', 'Inches', 'Position', [10 2 3 3])

for e = 154
    clf;
    f = 5;
    a = 4;

    axes('Position', [.2 .55 .75 .4]); hold on
    t = FoldedPSTH(spikeCounts{f, a, e}, binEdges, 5, 5);
    AlphaLine(binCenter, t, motorColor)
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', [])
    yLims = get(gca, 'YLim');
%     text(-.9, yLims(2), ColorText('Area 1', humanColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    text(-.9, yLims(2), ColorText('Motor', motorColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

    axes('Position', [.2 .1 .75 .4]); hold on
    r = Rasterfy(spikeCounts{f, a, e}, .4, true);
    plot(r(:,1), r(:,2), 'Color', motorColor, 'LineWidth', 1)
    set(gca, 'YLim', [0 length(spikeCounts{f,a,e})], 'YTick', [], 'YColor', 'none')
    
    xlabel('Time (seconds)')
    print(gcf, fullfile('.\Figures\SFN', sprintf('HumanIndent_Motor%d', e)), '-dpng', '-r900')
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
% humanColor = rgb(67, 160, 71);

motorColor = rgb(0, 137, 123);
% close all;
figure('Units', 'Inches', 'Position', [8 4 3 2])
for e = 161
    clf;    
    axes('Position', [.175 .25 .8 .715]); hold on
    for f = 3 
            t = FoldedPSTH(vertcat(spikeCounts{f, :, e}), binEdges, 5, 5);
            AlphaLine(binCenter, t, motorColor)
    end
    ylabel('Firing Rate (Hz)')
    set(gca, 'XTick', -1:1:2)
    xlabel('Time (s)')
    
    yLims = get(gca, 'YLim');
    text(-.9, yLims(2), ColorText('Motor', motorColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')


%     % Stimulus lines
%     plot([0 0], [25 40], '--', 'Color', [200, 200, 200] ./ 255)
%     plot([1 1], [25 40], '--', 'Color', [200, 200, 200] ./ 255)

    print(gcf, fullfile('.\Figures\SFN', sprintf('HumanIndent_ME%d', e)), '-dpng', '-r900')
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
for e = 161
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




%% Load Cortical macaque data

load('S:\UserFolders\CharlesGreenspon\SpeciesVibrations\RawMacaqueCorticalIndentationData.mat')

a1Idx = find(strcmp(L4_str.specs.anat_loc, 'a1'));
macaqueAI_area1 = L4_str.data.transienceIndex(a1Idx);

a3bIdx = find(strcmp(L4_str.specs.anat_loc, '3b'));
macaqueAI_area3b = L4_str.data.transienceIndex(a3bIdx);

nhpColor = rgb(251, 140, 0);

%%
% Example Rasters - Macaque
SetFont('Arial', 14)

dt = 0.02;
binEdges = 0:dt:1;
binCenter = binEdges(1:end-1) + (dt / 2);


% figure('Units', 'Inches', 'Position', [10 2 3 3])

% coi = [14, 72, 117];
for e = a1Idx'
    clf;

    axes('Position', [.2 .55 .75 .4]); hold on
    t = FoldedPSTH(cellfun(@transpose, L4_str.data.rasters{e}, 'UniformOutput', false), binEdges, 5, 5);
    AlphaLine(binCenter, t, nhpColor)
    ylabel('Firing Rate (Hz)')
    set(gca, 'XLim', [0 1], 'XTick', [])
    yLims = get(gca, 'YLim');

%     if strcmpi(L4_str.specs.anat_loc{e}, '3B')
%         text(1, yLims(2), ColorText('Area 3b', nhpColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
%     else
        text(1, yLims(2), ColorText('Area 1', nhpColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
%     end

    axes('Position', [.2 .1 .75 .4]); hold on
    r = Rasterfy(cellfun(@transpose, L4_str.data.rasters{e}, 'UniformOutput', false), .4, true);
    plot(r(:,1), r(:,2), 'Color', nhpColor)
    set(gca, 'XLim', [0 1], 'YTick', [], 'YColor', 'none')

    %     set(gca, 'XLim', [-1 2], 'YLim', [0 length(r(~isnan(r(:, 1))))], 'YTick', [], 'YColor', 'none')
%     xlabel('Time (seconds)')
    
    print(gcf, fullfile('.\Figures\Macaque', sprintf('Macaque_%s_%d', L4_str.specs.anat_loc{e}, e)), '-dpng', '-r900')

end

%% Only psth

% figure('Units', 'Inches', 'Position', [8 4 3 2])

dt = 0.02;
binEdges = 0:dt:1;
binCenter = binEdges(1:end-1) + (dt / 2);

for e = 107
    clf;    
    axes('Position', [.175 .25 .8 .715]); hold on
       t = FoldedPSTH(cellfun(@transpose, L4_str.data.rasters{e}, 'UniformOutput', false), binEdges, 5, 5);
    AlphaLine(binCenter, t, nhpColor)
    ylabel('Firing Rate (Hz)')
    set(gca, 'XLim', [0 1])
    yLims = get(gca, 'YLim');

    if strcmpi(L4_str.specs.anat_loc{e}, '3B')
        text(1, yLims(2), ColorText('Area 3b', nhpColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
    else
        text(1, yLims(2), ColorText('Area 1', nhpColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
    end

    
    yLims = get(gca, 'YLim');
    text(-.9, yLims(2), ColorText('Motor', motorColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

    print(gcf, fullfile('.\Figures\SFN', sprintf('Macaque_%s_%d', L4_str.specs.anat_loc{e}, e)), '-dpng', '-r900')
end





% figure('Units', 'Inches', 'Position', [10 2 4 6])
% r = Rasterfy(indentation_data.aligned_spikes(1, :), .4, true);
% plot(r(:,1), r(:,2), 'Color', nhpColor)

%% Load Cuneate macaque data






%% Adaptibility Index Figure

nhpColor = rgb(251, 140, 0);
humanColor = rgb(57, 73, 171);
motorColor = rgb(0, 137, 123);

sensoryChIdx = [65:96, 193:224]; sigSensory = sensoryChIdx(ismember(sensoryChIdx, sigElecIdx));
motorChIdx = [1:64, 97:192, 225:256]; sigMotor = motorChIdx(ismember(motorChIdx, sigElecIdx));

% figure('Units', 'Inches', 'Position', [10 2 4 8])
figure('Units', 'Inches', 'Position', [8 4 3 2.5])

clf;
% axes('Position', [.15 .8 .8 .2])
histogram(macaqueAI_area3b, 'BinEdges', linspace(0, 1, 7), 'Normalization', 'probability', ...
    'FaceColor', nhpColor)
box off
set(gca, 'YLim', [0 .3], 'YTick', 0:.1:.6, 'XTick', [])
text(0, .3, ColorText({'Macaque - Area 3b'}, nhpColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
xlabel('Adaptation Index')
% ylabel('Proportion')
print(gcf, fullfile('.\Figures\SFN', 'Macaque_3B_adaptationIndex'), '-dpng', '-r900')


axes('Position', [.15 .55 .8 .2])
histogram(macaqueAI_area1, 'BinEdges', linspace(0, 1, 7), 'Normalization', 'probability', ...
    'FaceColor', nhpColor)
box off
set(gca, 'YLim', [0 .6], 'YTick', 0:.2:.6, 'XTick', [])
text(0, .3, ColorText({'Macaque - Area 1'}, nhpColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
print(gcf, fullfile('.\Figures\SFN', 'Macaque_A1_adaptationIndex'), '-dpng', '-r900')


axes('Position', [.15 .3 .8 .2])
histogram(adaptationIndex(sigSensory), 'BinEdges', linspace(0, 1, 7), 'Normalization', 'probability', ...
    'FaceColor', humanColor)
box off
set(gca, 'YLim', [0 .6], 'YTick', 0:.2:.6, 'XTick', [])
text(0, .6, ColorText({'Human - Area 1'}, humanColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
print(gcf, fullfile('.\Figures\SFN', 'Human_A1_adaptationIndex'), '-dpng', '-r900')



axes('Position', [.15 .05 .8 .2])
histogram(adaptationIndex(sigMotor), 'BinEdges', linspace(0, 1, 7), 'Normalization', 'probability', ...
    'FaceColor', motorColor)
box off
set(gca, 'YLim', [0 .6], 'YTick', 0:.2:.6, 'XTick', [])
text(0, .3, ColorText({'Human - Motor'}, motorColor), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
print(gcf, fullfile('.\Figures\SFN', 'Human_MC_adaptationIndex'), '-dpng', '-r900')


xlabel('Adaptation Index')
ylabel('Proportion', 'Position', [-.2 .325 -1])
   
print(gcf, fullfile('.\Figures\SFN', 'adaptationIndex'), '-dpng', '-r900')
