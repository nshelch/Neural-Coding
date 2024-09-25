function [phase_lock_ratio, isi_phase_distribution] = ISIPhaseLock(spike_times, frequency, discretize)
    if nargin == 2
        discretize = true;
    end
    % Ensure inputs are row vectors
    spike_times = cellfun(@(c) c(:), spike_times, 'UniformOutput', false);
    % Get inter-spike-interval for all combinations of spikes for each trial
    isis = cellfun(@(c) pdist(c), spike_times, 'UniformOutput', false);
    % Concatenate across trials
    isis = cat(2, isis{:});
    % Remove 0 distance isis (artifact from using pdist to compute isi)
    isis(isis == 0) = [];
    % Normalize the isis by the reciprocal of the frequency
    isi_phase_distribution = isis / (1/frequency);
    % Find offset relative to phase (round lets you compute phase across multiple stimulus phases)
    isi_phase_distribution = isi_phase_distribution - round(isi_phase_distribution);
    % Compute ratio of high phase locking (within 1/6 pi) isis
    phase_lock_ratio = sum(isi_phase_distribution < 1/12 & isi_phase_distribution > -1/12) / length(isi_phase_distribution);
    if discretize
        isi_phase_distribution = histcounts(isi_phase_distribution, linspace(-0.5, .5, 12), 'Normalization', 'probability');
    end
end