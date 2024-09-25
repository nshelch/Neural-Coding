load("U:\UserFolders\CharlesGreenspon\SpeciesVibrations\RawSpeciesVibData.mat")
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

%% First need to filter responsive neurons
for i = 1:height(all_data)
    [usp, ~, u_idx] = unique(all_data{i, "StimParams"}{1}, "rows");

    for p = 1:height(usp)
        p_idx = u_idx == p;
        temp_spikes = {all_data{i, "SpikeTimes"}{1}{p_idx}};
        temp_spikes = cellfun(@(c) c(c > 0), temp_spikes, 'UniformOutput', false);
        spike_counts = cellfun();
    end
    all_data{i, "SpikeRates"} = {usp};
end


%% Rate intensity functions



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
ur = unique(all_data(:, ["Species", "Location"]), 'rows');
% Sort by heirarchy
% idx = [{5,4,3}, {2}, {7,6}, {1}]; % Cells allow multiple on a given plot
clf; hold on
idx = find(strcmpi(all_data.Species, 'Macaque') & strcmpi(all_data.Location, 'SA'));
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
    print(gcf, fullfile(export_path, sprintf('%s_%s_%d', all_data.Species{i}, all_data.Location{i}, i)),...
        '-dpng', '-r300')
end