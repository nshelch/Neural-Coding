SetFont('Arial', 15)

%% Vibration Rasters + Rate Intensity Functions

% 59: Human area 1 - ch 219
% 91: Macaque area 1 - ch 27
% 101: Macaque area 3b - ch 37
% 195: Macaque cuneate - ch 49
coi = [59, 91, 101, 195];

figure('Units', 'Inches', 'Position', [10 2 4 8])

all_freq = [unique(all_data{coi(1), 'RateIntensity'}{1}{:, "Freq"}); unique(all_data{coi(2), 'RateIntensity'}{1}{:, "Freq"}); ...
    unique(all_data{coi(3), 'RateIntensity'}{1}{:, "Freq"}); unique(all_data{coi(4), 'RateIntensity'}{1}{:, "Freq"})];

uniq_freq = unique(all_freq);
freqColor = cmap_gradient([rgb(3, 155, 229); rgb(216, 27, 96)], length(uniq_freq));

for i = 59

    % Raster
    [usp, ~, u_idx] = unique(all_data{i, "StimParams"}{1}, "rows");

    clf;
    axes('Position', [.175 .375 .775 .575]); hold on
    r0 = 0;
    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{i, "SpikeTimes"}{1}{p_idx}};
        r = Rasterfy(temp_spikes, .4, true);

        % Set the color based on the frequency
        fcIdx = uniq_freq == usp{p, "Freq"} ;
        plot(r(:,1), r(:,2) + r0, 'Color', freqColor(fcIdx,:))
        r0 = r0 + length(temp_spikes);
    end
    xlabel('Time (s)');
    set(gca, 'YLim', [0 r0], ...
        'XLim', [-.5 1.5], ...
        'YTick', [], 'YColor', 'none')

    freq = unique(all_data{i, 'RateIntensity'}{1}{:, "Freq"});
    for f = 1:length(freq)
        fcIdx = freq(f) == uniq_freq;
        %         text(-.25, 50*f - 10, ColorText(sprintf('%i Hz', freq(f)), freqColor(fcIdx,:)), 'HorizontalAlignment', 'right') % Macaque
        text(-.5, 150*f - 50, ColorText(sprintf('%i Hz', freq(f)), freqColor(fcIdx,:)), 'HorizontalAlignment', 'right')
    end

    title(sprintf('%s: %s', all_data.Species{i}, all_data.Location{i}))

    % Rate Intensity
    axes('Position', [.175 .085 .775 .2]); hold on
    for f = 1:length(freq)
        fcIdx = freq(f) == uniq_freq;

        freqIdx = all_data{i, 'RateIntensity'}{1}{:, "Freq"} == freq(f);
        temp_data = vertcat(all_data{i, 'RateIntensity'}{1}{freqIdx, "SpikeCount"}{:})';

        AlphaLine(all_data{i, 'RateIntensity'}{1}{freqIdx, "Amp"}, temp_data, freqColor(fcIdx, :), 'ErrorType', 'SEM')
        hold on
    end

    set(gca, 'XScale', 'log', 'XLim', [0 1000], 'XTick', [0 10 100 1000])
    ylabel('Firing Rate [Hz]')
    xlabel('Amplitude')

    print(gcf, fullfile('.\Figures\SFN', sprintf('Vibrations_%s_%s_Ch%d', all_data.Species{i}, all_data.Location{i}, all_data.CellID(i))), '-dpng', '-r900')

end

%% Tuning Curves
coi = [59, 91, 101, 195];

figure('Units', 'Inches', 'Position', [10 2 4 3])
nhpColor = rgb(251, 140, 0);
humanColor = rgb(57, 73, 171); %rgb(67, 160, 71);

for i = [59, 91, 101, 195]

    x = unique(all_data.RateIntensity{i}.Freq);
    y = all_data.freqResponse{i};  % Replace with your actual response data

    gammaModel = fittype('A * x^(k-1) * exp(-x/theta) + C', ...
        'independent', 'x', ...
        'coefficients', {'A', 'k', 'theta', 'C'});

    %initial guesses for A, k, and theta
    startPoints = [1, 2, 100, 0];

    %set lower and upper bounds for the parameters to avoid problematic values
    lowerBounds = [0, 0.1, 1, -Inf];  % Lower bounds for A, k, theta, C
    upperBounds = [Inf, Inf, Inf, Inf];  % Upper bounds for A, k, theta, C

    %fit
    gammaFit = fit(x, y, gammaModel, 'StartPoint', startPoints, 'Lower', lowerBounds, 'Upper', upperBounds);

    %evaluate the gamma fit over the range of x values
    x_fit = linspace(min(x), max(x), 1000);

    %evaluate the gamma fit at these x values
    y_fit = feval(gammaFit, x_fit);

    %plot the data points in light grey
    clf;
    axes('Position', [.175 .2 .775 .7]); hold on

    if strcmpi(all_data.Species{i}, 'Human')
        plot(x, y, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', humanColor, 'MarkerFaceColor', humanColor);
    else
        plot(x, y, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', nhpColor, 'MarkerFaceColor', nhpColor);
    end

    hold on

    %plot the fitted gamma line in dark lilac
    plot(x_fit, y_fit, 'Color', 'k', 'LineWidth', 1.5);

    box off
    %set x-ticks to only the frequency values in x
    set(gca, 'XLim', [0 800], 'XTick', 0:200:800);

    title(sprintf('%s: %s', all_data.Species{i}, all_data.Location{i}))
    xlabel('Frequency [Hz]');
    ylabel('Firing Rate [Hz]');
    print(gcf, fullfile('.\Figures\SFN', sprintf('TuningCurve_%s_%s_Ch%d', all_data.Species{i}, all_data.Location{i}, all_data.CellID(i))), '-dpng', '-r900')

end

%% Rate Heatmaps per area

%% Tuning curves across areas

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'SA');
sa_TC = all_data{idx, "whm"};

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'RA');
ra_TC = all_data{idx, "whm"};

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'PC');
pc_TC = all_data{idx, "whm"};

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Cuneate');
cun_TC = all_data{idx, "whm"};

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Area 3b');
area3b_TC = all_data{idx, "whm"};

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Area 1');
area1_TC = all_data{idx, "whm"};

idx = strcmpi(all_data.Species, 'Human') & strcmpi(all_data.Location, 'Area 1');
human_TC = all_data{idx, "whm"};

% Figure
nhpColor = rgb(251, 140, 0);
humanColor = rgb(57, 73, 171); %rgb(67, 160, 71);
markerSize = 10;
xPos = [1,2,3, 5, 7,8, 10];

% figure('Units', 'Inches', 'Position', [10 2 6 4])
clf; hold on

Swarm(xPos(1), sa_TC, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50)
Swarm(xPos(2), ra_TC, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50)
Swarm(xPos(3), pc_TC, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50)
Swarm(xPos(4), cun_TC, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)
Swarm(xPos(5), area3b_TC, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)
Swarm(xPos(6), area1_TC, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)
Swarm(xPos(7), human_TC, 'Color', humanColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)

set(gca, 'YLim', [-10 600], 'XLim', [0.5 11], 'XTick', xPos, 'XTickLabel', {'SA', 'RA', 'PC', 'Cuneate', '3b', '1', '1'}, ...
    'XTickLabelRotation', 0)
ylabel('Tuning Width')

text(0.75, 600, ColorText({'Macaque', 'Human'}, [nhpColor; humanColor]), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

print(gcf, fullfile('.\Figures\SFN', 'TuningAcrossAreas'), '-dpng', '-r900')


%% Macaque phase lock across areas

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'SA');
sa_PL = cellfun(@(z) nanmean(z{:, "PhaseLockRatio"}), all_data{idx, "PhaseLocking"});

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'RA');
ra_PL = cellfun(@(z) nanmean(z{:, "PhaseLockRatio"}), all_data{idx, "PhaseLocking"});

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'PC');
pc_PL = cellfun(@(z) nanmean(z{:, "PhaseLockRatio"}), all_data{idx, "PhaseLocking"});

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Cuneate');
cun_PL = cellfun(@(z) nanmean(z{:, "PhaseLockRatio"}), all_data{idx, "PhaseLocking"});

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Area 3b');
area3b_PL = cellfun(@(z) nanmean(z{:, "PhaseLockRatio"}), all_data{idx, "PhaseLocking"});

idx = strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'Area 1');
area1_PL = cellfun(@(z) nanmean(z{:, "PhaseLockRatio"}), all_data{idx, "PhaseLocking"});

idx = strcmpi(all_data.Species, 'Human') & strcmpi(all_data.Location, 'Area 1');
human_PL = cellfun(@(z) nanmean(z{:, "PhaseLockRatio"}), all_data{idx, "PhaseLocking"});

% Figure
nhpColor = rgb(251, 140, 0);
humanColor = rgb(57, 73, 171); %rgb(67, 160, 71);
markerSize = 10;
xPos = [1,2,3, 5, 7,8, 10];

figure('Units', 'Inches', 'Position', [10 2 6 4])
clf; hold on
% scatter(1, sa_PL, 20, nhpColor, 'filled')
% scatter(1.5, ra_PL, 20, nhpColor, 'filled')
% scatter(2, pc_PL, 20, nhpColor, 'filled')
% scatter(3, cun_PL, 20, nhpColor, 'filled')
% scatter(4, area3b_PL, 20, nhpColor, 'filled')
% scatter(4.5, area1_PL, 20, nhpColor, 'filled')
% scatter(5.5, human_PL, 20, humanColor, 'filled')

Swarm(xPos(1), sa_PL, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50)
Swarm(xPos(2), ra_PL, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50)
Swarm(xPos(3), pc_PL, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50)
Swarm(xPos(4), cun_PL, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)
Swarm(xPos(5), area3b_PL, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)
Swarm(xPos(6), area1_PL, 'Color', nhpColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)
Swarm(xPos(7), human_PL, 'Color', humanColor, 'DS', 'Violin', 'CLW', 1.5, 'SMS', markerSize, 'SPL', 50, 'DW', .5)

set(gca, 'YLim', [0 1], 'YTick', 0:.2:1, 'XLim', [0.5 11], 'XTick', xPos, 'XTickLabel', {'SA', 'RA', 'PC', 'Cuneate', '3b', '1', '1'}, ...
    'XTickLabelRotation', 0)
ylabel('Phase Locking Probability')

text(0.75, 1, ColorText({'Macaque', 'Human'}, [nhpColor; humanColor]), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

print(gcf, fullfile('.\Figures\SFN', 'PhaseLockAcrossAreas'), '-dpng', '-r900')


%% Phase locking per area

freq = [5, 10, 20, 50, 75, 100, 150, 200, 250, 300, 400, 600, 800];

% Figure
nhpColor = rgb(251, 140, 0);
humanColor = rgb(57, 73, 171); %rgb(67, 160, 71);
markerSize = 5;

figure('Units', 'Inches', 'Position', [10 2 4 8])

clf;

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

axes('Position', [.175 .75 .775 .2]); hold on
scatter(freq, y, 20, nhpColor, 'filled')
text(800, 1, ColorText({'Cuneate'}, nhpColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

xlim([0 800])
ylim([0 1])
set(gca, 'XTickLabel', '')


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

axes('Position', [.175 .525 .775 .2]); hold on
scatter(freq, y, 20, nhpColor, 'filled')
ylim([0 1])
set(gca, 'XTickLabel', '')
text(800, 1, ColorText({'Area 3b'}, nhpColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')

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

axes('Position', [.175 .3 .775 .2]); hold on
scatter(freq, y, 20, nhpColor, 'filled')
ylim([0 1])
set(gca, 'XTickLabel', '')
text(800, 1, ColorText({'Area 1'}, nhpColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')


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

axes('Position', [.175 .075 .775 .2]); hold on
scatter(freq, y, 20, humanColor, 'filled')
ylim([0 1])
xlim([0 800])
text(800, 1, ColorText({'Area 1'}, humanColor), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
xlabel('Frequency [Hz]')
ylabel('Phase Locking Probability')

print(gcf, fullfile('.\Figures\SFN', 'PhaseLockWithinAreas'), '-dpng', '-r900')
