%% Human Area 1
load('U:\UserFolders\CharlesGreenspon\BCI_TactileCoding\ProcessedData\Vibrations_data')
% Filter only sensory electrodes
cmap = LoadSubjectChannelMap('BCI02');
ch_filter = cat(1, cmap.ArrayLocations{cmap.IsSensory});
ch_filter = sort(ch_filter(~isnan(ch_filter)));
num_neurons = length(ch_filter);
species = repmat({'Human'}, [num_neurons, 1]);
location = repmat({'area 1'}, [num_neurons, 1]);
all_spikes = cat(1, data.Spikes{:}); all_spikes = all_spikes(:, ch_filter);
[spike_times, stim_params] = deal(cell(num_neurons, 1));
for i = 1:num_neurons
    spike_times{i} = all_spikes(:,i);
    stim_params{i} = data(:, ["Freq", "Amp"]);
end
HumanCorticalSineData = table(species, location, ch_filter, spike_times, stim_params, ...
    'VariableNames', ["Species", "Location", "CellID", "SpikeTimes", "StimParams"]);

%% Macaque cortical sine
load('V:\ProjectFolders\ArchivedProjects\TerminalCuneate\Data\Sine\CorticalSineDataHigh_Raw.mat')
num_neurons = length(L2_str.neuron_id);
cell_id = [1:num_neurons]';
species = repmat({'Macaque'}, [num_neurons, 1]);
location = L2_str.anat_loc;
[spike_times, stim_params] = deal(cell(num_neurons, 1));
for i = 1:num_neurons
    spike_times{i} = L2_str.ObservedSpikes{i}';
    sp = L2_str.state_v{i}(:, [1,2]);
    stim_params{i} = table(sp(:,1), sp(:,2), ...
        'VariableNames', {'Freq', 'Amp'});
end
MacaqueCorticalSineData = table(species, location', cell_id, spike_times, stim_params, ...
    'VariableNames', ["Species", "Location", "CellID", "SpikeTimes", "StimParams"]);

%% Macaque peripheral Chubbeck
load('V:\ProjectFolders\ArchivedProjects\S1Texture\mat\pma_ChubSine1.mat')
num_neurons = length(L2_str.neuron_id);
cell_id = [1:num_neurons]';
species = repmat({'Macaque'}, [num_neurons, 1]);
location = L2_str.nType;
[spike_times, stim_params] = deal(cell(num_neurons, 1));
for i = 1:num_neurons
    spike_times{i} = L2_str.ObservedSpikes{i}(:);
    sp = repmat(L2_str.state_v{i}(:, [4,3]), [size(L2_str.ObservedSpikes{i}, 2), 1]);
    stim_params{i} = table(sp(:,1), sp(:,2), ...
        'VariableNames', {'Freq', 'Amp'});
end
MacaqueAfferentSineData = table(species, location, cell_id, spike_times, stim_params, ...
    'VariableNames', ["Species", "Location", "CellID", "SpikeTimes", "StimParams"]);

%% Macaque Cuneate UMEA
load('V:\ProjectFolders\ArchivedProjects\TerminalCuneate\Data\Sine\ProcessedSineData.mat')
num_neurons = size(sine_data.spike_times, 1);
species = repmat({'Macaque'}, [num_neurons, 1]);
location = repmat({'Cuneate'}, [num_neurons, 1]);
cell_id = [1:num_neurons]';

[spike_times, stim_params] = deal(cell(num_neurons, 1));
for i = 1:num_neurons
    spike_times{i} = sine_data.spike_times(i,:)';
    stim_params{i} = table(sine_data.stim_params{i}(:,1), sine_data.stim_params{i}(:,2), ...
        'VariableNames', {'Freq', 'Amp'});
end
MacaqueCuneateSineData = table(species, location, cell_id, spike_times, stim_params, ...
    'VariableNames', ["Species", "Location", "CellID", "SpikeTimes", "StimParams"]);

%% Save and export
clearvars -except *SineData
