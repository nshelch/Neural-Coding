%% Human Area 1
load('U:\UserFolders\CharlesGreenspon\BCI_TactileCoding\ProcessedData\Indentations_data')

% Filter only sensory electrodes
chanMap = LoadSubjectChannelMap('BCI02');

chFilter = cat(1, chanMap.ArrayLocations{chanMap.IsSensory});
chFilter = sort(chFilter(~isnan(chFilter)));
numNeurons = length(chFilter);

species = repmat({'Human'}, [numNeurons, 1]);
location = repmat({'area 1'}, [numNeurons, 1]);

allSpikes = cat(1, data.Spikes{:}); allSpikes = allSpikes(:, chFilter);
[spikeTimes, stimParams] = deal(cell(numNeurons, 1));
for i = 1:numNeurons
    spikeTimes{i} = allSpikes(:,i);
    stimParams{i} = data(:, ["Freq", "Amp"]);
end

HumanCorticalIndData = table(species, location, chFilter, spikeTimes, stimParams, ...
    'VariableNames', ["Species", "Location", "CellID", "SpikeTimes", "StimParams"]);

%% Macaque cortical sine
load('U:\UserFolders\CharlesGreenspon\SpeciesVibrations\RawMacaqueCorticalIndentationData.mat')
numNeurons = length(L4_str.specs.neuron_id);
cellID = [1:numNeurons]';
species = repmat({'Macaque'}, [numNeurons, 1]);
location = L4_str.specs.anat_loc;

[spikeTimes, stimParams] = deal(cell(numNeurons, 1));
for i = 1:numNeurons
    spikeTimes{i} = L2_str.ObservedSpikes{i}';
    sp = L2_str.state_v{i}(:, [1,2]);
    stimParams{i} = table(sp(:,1), sp(:,2), ...
        'VariableNames', {'Freq', 'Amp'});
end
MacaqueCorticalIndData = table(species, location', cellID, spikeTimes, stimParams, ...
    'VariableNames', ["Species", "Location", "CellID", "SpikeTimes", "StimParams"]);
