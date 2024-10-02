load("S:\UserFolders\CharlesGreenspon\SpeciesVibrations\RawSpeciesVibData.mat")
addpath(genpath('HelperFunctions'))
all_data = cat(1, HumanCorticalSineData, MacaqueCorticalSineData, MacaqueCuneateSineData, MacaqueAfferentSineData);
clearvars -except all_data
warning('off', 'MATLAB:table:RowsAddedNewVars')

%% Example rasters for each
id = [60, 79, 77, 149, 229, 230]; %68
fcol = cool();
freq_idx = linspace(log(10), log(500), 256);
clf; 
tiledlayout('flow');
for i = 1:length(id)
    nexttile(); hold on
    r0 = 0;
    [usp, ~, u_idx] = unique(all_data{id(i), "StimParams"}{1}, "rows");
    uf = unique(usp.Freq);
    
    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{id(i), "SpikeTimes"}{1}{p_idx}};
        r = Rasterfy(temp_spikes, .4, true);
        [~,c_idx] = min(abs(freq_idx - log(usp.Freq(p))));
        plot(r(:,1), r(:,2) + r0, 'Color', fcol(c_idx,:))
        r0 = r0 + length(temp_spikes);
    end

    xlabel('Time (s');
    set(gca, 'YLim', [0 r0], ...
             'XLim', [-.5 1.5], ...
             'YTick', [])
    title(sprintf('%s: %s', all_data.Species{id(i)}, all_data.Location{id(i)}))

end

%% First need to filter responsive neurons (FINISH)

baseline_int = [-.5 -.3];
baseline2_int = [-.3 -.1];
stim_int = [0 1];

for i = 1:height(all_data)
    [usp, ~, u_idx] = unique(all_data{i, "StimParams"}{1}, "rows");

    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{i, "SpikeTimes"}{1}{p_idx}};
%         temp_spikes = cellfun(@(c) c(c > 0), temp_spikes, 'UniformOutput', false);
        baseline_spikes = cellfun(@(c) sum(c >= baseline_int(1) & c < baseline_int(2)), temp_spikes) ./ diff(baseline_int);
        baseline2_spikes = cellfun(@(c) sum(c >= baseline2_int(1) & c < baseline2_int(2)), temp_spikes) ./ diff(baseline2_int);
        stim_spikes = cellfun(@(c) sum(c >= stim_int(1) & c <= stim_int(1)), temp_spikes) ./ diff(stim_int);
        usp{p, "BaselineFR"} = {baseline_spikes};
        usp{p, "Baseline2FR"} = {baseline2_spikes};
        usp{p, "StimulusFR"} = {stim_spikes};
    end

    % Permutation


    all_data{i, "SpikeRates"} = {usp};
end


%% Rate intensity functions
for i = 1:64 %height(all_data)
    [usp, ~, u_idx] = unique(all_data{i, "StimParams"}{1}, "rows");
    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{i, "SpikeTimes"}{1}{p_idx}};
        stim_spikes = cellfun(@(c) sum(c >= 0 & c <= 1), temp_spikes); % Maybe extend it to 1.2 for human data
        usp{p, "SpikeCount"} = {stim_spikes};
        usp{p, "AvgFR"} = nanmean(stim_spikes);
    end
    all_data{i, "RateIntensity"} = {usp};

    % Get the best response across amplitudes
    freq = unique(usp{:, "Freq"});
    tmp_fr = zeros(size(freq));
    for f = 1:length(freq)
        f_idx = usp{:, "Freq"} == freq(f);
        tmp_fr(f) = prctile(usp{f_idx, "AvgFR"}, 75);
    end

    plot(freq, tmp_fr, 'k')
    % Tuning curves for rate modulation stuff

    p = polyfit(freq, tmp_fr, 6); % Lee 24 paper used 6 for their polyfit
% 
    y1 = polyval(p, linspace(min(freq), max(freq), 100));
    hold on
    plot(linspace(min(freq), max(freq), 100), y1)
    hold off
    %
    %     usf = table(freq, 'VariableNames', {'Freq'});

        print(gcf, fullfile('.\Figures\Vibrations\TuningCurves', sprintf('%s_%s_Ch%d', all_data.Species{i}, all_data.Location{i}, all_data.CellID(i))), '-dpng', '-r900')


end

%% Figure plotting
figure('Units', 'Inches', 'Position', [10 2 4 8])

fcol = cool();
freq_idx = linspace(log(10), log(500), 256);

for i = 65:height(all_data)

    freq = unique(all_data{i, 'RateIntensity'}{1}{:, "Freq"});
    amp = unique(all_data{i, 'RateIntensity'}{1}{:, "Amp"});
    freqColor = cmap_gradient([rgb(3, 155, 229); rgb(216, 27, 96)], length(freq));

    % Raster
    r0 = 0;
    [usp, ~, u_idx] = unique(all_data{i, "StimParams"}{1}, "rows");

    clf;
    axes('Position', [.15 .375 .8 .6]); hold on
    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{i, "SpikeTimes"}{1}{p_idx}};
        r = Rasterfy(temp_spikes, .4, true);
        [~,c_idx] = min(abs(freq_idx - log(usp.Freq(p))));
        plot(r(:,1), r(:,2) + r0, 'Color', fcol(c_idx,:))
        r0 = r0 + length(temp_spikes);
    end
    xlabel('Time (s');
    set(gca, 'YLim', [0 r0], ...
        'XLim', [-.5 1.5], ...
        'YTick', [])
    title(sprintf('%s: Ch: %i', all_data.Species{i}, all_data.CellID(i)))

    % Rate Intensity
    axes('Position', [.15 .075 .8 .225]); hold on
    for f = 1:length(freq)
        freqIdx = all_data{i, 'RateIntensity'}{1}{:, "Freq"} == freq(f);
        temp_data = vertcat(all_data{i, 'RateIntensity'}{1}{freqIdx, "SpikeCount"}{:})';

        AlphaLine([5, 50, 100, 200, 400], temp_data, freqColor(f, :), 'ErrorType', 'SEM')
        hold on

    end
    set(gca, 'XScale', 'log', 'XLim', [0 400])
    ylabel('Firing Rate [Hz]')
    xlabel('Amplitude')

    print(gcf, fullfile('.\Figures\Vibrations\RateIntensity', sprintf('%s_%s_Ch%d', all_data.Species{i}, all_data.Location{i}, all_data.CellID(i))), '-dpng', '-r900')

end
%% Compute phase locking
for i = 1:height(all_data)
    [usp, ~, u_idx] = unique(all_data{i, "StimParams"}{1}, "rows");
    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{i, "SpikeTimes"}{1}{p_idx}};
        temp_spikes = cellfun(@(c) c(c > 0), temp_spikes, 'UniformOutput', false);
        [phase_lock_ratio, isi_phase_distribution] = ISIPhaseLock(temp_spikes, usp{p, "Freq"});
        usp{p, "PhaseLockRatio"} = phase_lock_ratio;
        usp{p, "ISIPhaseDist"} = {isi_phase_distribution};
    end
    all_data{i, "PhaseLocking"} = {usp};
end

%% Plot

%% Macaque phase lock across areas
nhpColor = rgb(251, 140, 0);

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'SA'));
sa_PL = zeros(size(idx));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    sa_PL(i) = nanmean(pls{:, "PhaseLockRatio"});
end

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'RA'));
ra_PL = zeros(size(idx));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    ra_PL(i) = nanmean(pls{:, "PhaseLockRatio"});
end

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'PC'));
pc_PL = zeros(size(idx));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    pc_PL(i) = nanmean(pls{:, "PhaseLockRatio"});
end

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Cuneate'));
cun_PL = zeros(size(idx));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    cun_PL(i) = nanmean(pls{:, "PhaseLockRatio"});
end

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Area 3b'));
area3b_PL = zeros(size(idx));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    area3b_PL(i) = nanmean(pls{:, "PhaseLockRatio"});
end

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Area 1'));
area1_PL = zeros(size(idx));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    area1_PL(i) = nanmean(pls{:, "PhaseLockRatio"});
end

clf; hold on
scatter(1, sa_PL, 20, nhpColor, 'filled')
scatter(2, ra_PL, 20, nhpColor, 'filled')
scatter(3, pc_PL, 20, nhpColor, 'filled')
scatter(4, cun_PL, 20, nhpColor, 'filled')
scatter(5, area3b_PL, 20, nhpColor, 'filled')
scatter(6, area1_PL, 20, nhpColor, 'filled')
set(gca, 'XLim', [0.5 6.5], 'XTick', 1:6, 'XTickLabel', {'SA', 'RA', 'PC', 'Cuneate', 'Area 3b', 'Area 1'})
ylabel('Phase locking probability')

idx = find(strcmpi(all_data.Species, 'Human') & strcmpi(all_data.Location, 'Area 1'));
human_PL = zeros(size(idx));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    human_PL(i) = nanmean(pls{:, "PhaseLockRatio"});
end

humanColor = rgb(67, 160, 71);
scatter(7, human_PL, 20, humanColor, 'filled')
set(gca, 'XLim', [0.5 7.5], 'XTick', 1:7, 'XTickLabel', {'SA', 'RA', 'PC', 'Cuneate', 'Area 3b', 'Area 1', 'Area 1 (H)'})

text(7.5, .8, ColorText({'Macaque', 'Human'}, [nhpColor; humanColor]), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

%% Macaque, Phase-lock per area
clf; 
freq = [5, 10, 20, 50, 75, 100, 150, 200, 250, 300, 400, 600, 800];

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Cuneate'));
y = NaN(length(idx), length(freq));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    for xi = 1:length(freq)
        xidx = pls.Freq == freq(xi);
        if sum(xidx) > 0
            y(i, xi) = nanmean(pls{xidx, "PhaseLockRatio"});
        end
    end
end
subplot(4, 1, 1)
scatter(freq, y, 20, nhpColor, 'filled')
title('Cuneate')
xlim([0 800])

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Area 3b'));
y = NaN(length(idx), length(freq));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    for xi = 1:length(freq)
        xidx = pls.Freq == freq(xi);
        if sum(xidx) > 0
            y(i, xi) = nanmean(pls{xidx, "PhaseLockRatio"});
        end
    end
end
subplot(4, 1, 2)
scatter(freq, y, 20, nhpColor, 'filled')
title('Area 3b')
ylim([0 1])

idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'AREA 1'));
y = NaN(length(idx), length(freq));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    for xi = 1:length(freq)
        xidx = pls.Freq == freq(xi);
        if sum(xidx) > 0
            y(i, xi) = nanmean(pls{xidx, "PhaseLockRatio"});
        end
    end
end
subplot(4, 1, 3)
scatter(freq, y, 20, nhpColor, 'filled')
title('Area 1')

idx = find(strcmpi(all_data.Species, 'Human') & strcmpi(all_data.Location, 'AREA 1'));
y = NaN(length(idx), length(freq));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    for xi = 1:length(freq)
        xidx = pls.Freq == freq(xi);
        if sum(xidx) > 0
            y(i, xi) = nanmean(pls{xidx, "PhaseLockRatio"});
        end
    end
end
subplot(4, 1, 4)
scatter(freq, y, 20, humanColor, 'filled')
title('Area 1')

%%
% Phase locking distribution
histogram('BinEdges', linspace(-0.5, .5, 12), 'BinCounts', pls{1, "ISIPhaseDist"}{:})



% ur = unique(all_data(:, ["Species", "Location"]), 'rows');
% % Sort by heirarchy
% % idx = [{5,4,3}, {2}, {7,6}, {1}]; % Cells allow multiple on a given plot
% clf; hold on
% idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'SA'));
for i = 1:length(idx)
    pls = all_data{idx(i), "PhaseLocking"}{1};
    x = unique(pls.Freq);
    y = zeros(size(x));
    for xi = 1:length(x)
        xidx = pls.Freq == x(xi);
        y(xi) = max(pls{xidx, "PhaseLockRatio"});
    end
    plot(x,y, 'Color', [.6 .6 .6])
end

%% Rate modulation





%% Scatter plot of firing rate modulation versus phase locking


return
%% Image for every unit
export_path = 'C:\Users\somlab\Desktop\ExportedImages\Vibrations\FrequencyRaster';
fcol = cool();
freq_idx = linspace(log(10), log(500), 256);
for i = 1:height(all_data)
    clf; hold on
    r0 = 0;
    [usp, ~, u_idx] = unique(all_data{i, "StimParams"}{1}, "rows");
    uf = unique(usp.Freq);
    
    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{i, "SpikeTimes"}{1}{p_idx}};
        r = Rasterfy(temp_spikes, .4, true);
        [~,c_idx] = min(abs(freq_idx - log(usp.Freq(p))));
        plot(r(:,1), r(:,2) + r0, 'Color', fcol(c_idx,:))
        r0 = r0 + length(temp_spikes);
    end

    xlabel('Time (s');
    set(gca, 'YLim', [0 r0], ...
             'XLim', [-.5 1.5], ...
             'YTick', [])
    title(sprintf('%s: %s', all_data.Species{i}, all_data.Location{i}))
%     print(gcf, fullfile(export_path, sprintf('%s_%s_%d', all_data.Species{i}, all_data.Location{i}, i)),...
%         '-dpng', '-r300')
end