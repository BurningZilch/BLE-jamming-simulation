rng(2);

% Initialize constants
k = 1.38e-23;  % Boltzmann constant in J/K
T = 290;       % Temperature in Kelvin
snr_dB = 20;
P_signal = 0.001;    % Adjusted Signal power in watts
B = 2e6;            % Bandwidth in Hz
R_b = 125e3;        % Bit rate in bits per second
P_dBm = 10 * log10(P_signal / 0.001);
N0 = k * T;         % Noise density in W/Hz
N0_dBmHz = 10 * log10(N0) + 30;  % Noise density in dBm/Hz
P_noise = N0 * B;

fc = 2.426e6;       % Channel frequency
c = 3e8;            % Speed of light (m/s)
lambda = c / fc;    % Wavelength (m)
sps = 2;
sampleRate = sps * 125e3;
dataLength = 128;   % Includes header, payload, and CRC
phyMode = "BR";
simMode = "LE1M";
packetType = "DM1";
channelModel = "Raytracing Channel";
enableEqualizer = false;
bitsPerByte = 8;

% Generate node and locator positions
nodeDistance = 2; 
xRange = [-5.6, 5.6]; 
yRange = [-5.6, 5.6]; 
zRange = [0, 4]; 
numTestPoints = 20;

Posnodes = [...
    xRange(1) + (xRange(2) - xRange(1)) * rand(1, numTestPoints);
    yRange(1) + (yRange(2) - yRange(1)) * rand(1, numTestPoints);
    zRange(1) + (zRange(2) - zRange(1)) * rand(1, numTestPoints)];

Poslocators = [
    2.8, 3.0, -1.8, 5.0, -2.2, 1.4, 4.0, -2.8, 3.6, -1.2;
   -2.0, 4.2,  6.0, -1.6, -4.4, 3.0, 0.0, -3.4, 4.6,  2.4;
    4.0, 4.0, 2.4, 1.2, 3.6, -0.8, 6.6, -4.0, -1.0, 5.4];

Posjammer = [0; 0; 1]; % Jammer at (0, 0, 1) meters
jammerPower = 0; % dBm
jammerPower_W = 10^((jammerPower - 30)/10); % Convert dBm to Watts


% Create scenario with generated positions
[locators, nodes, antPosNode] = createScenario(Posnodes, Poslocators);
[jammers, antPosJammer] = createJammer(Posjammer);

% Initialize impairments
initImpTemplate = helperBLEImpairmentsInit("LE Coded PHY", sps);

% Create PHY layers for locators and nodes
numLocators = size(Poslocators, 2);
numNodes = size(Posnodes, 2);
txPhys = cell(1, numLocators);
rxPhys = cell(1, numNodes);

for i = 1:numLocators
    txPhy = helperBluetoothPHY('TxPower', P_dBm, 'Mode', phyMode);
    txPhy.NodePosition = locators(i).AntennaPosition;
    txPhys{i} = txPhy;
end

for i = 1:numNodes
    rxPhy = helperBluetoothPHY('Mode', phyMode);
    rxPhy.NodePosition = nodes(i).AntennaPosition;
    rxPhys{i} = rxPhy;
end

% Initialize channels between locators and nodes
for i = 1:numNodes
    for a = 1:numLocators
        BLEChannel = helperBluetoothChannelInit(sampleRate, channelModel, locators(a).AntennaPosition, nodes(i).AntennaPosition);
        rxPhys{i}.BLEChannels{a} = BLEChannel;
    end
end

% Initialize jammer and its channels to nodes
jammerTx = helperBluetoothPHY('TxPower', jammerPower, 'Mode', phyMode);
jammerTx.NodePosition = Posjammer;



jammerChannels = cell(1, numNodes);
for i = 1:numNodes
    jammerChannels{i} = helperBluetoothChannelInit(sampleRate, 'Raytracing Channel', Posjammer, nodes(i).AntennaPosition);
end

viewer = siteviewer(SceneModel=BLEChannel.VisualVar.MapFileName);
% 
 show(nodes,Icon="bleRxIcon.png");
 show(locators,Icon="bleTxIcon.png");
 show(jammers,Icon="attacker.png");
rays_jammers = cell(1, length(jammers));  % Preallocate for jammers' rays
for j = 1:length(jammers)
    for i = 1:length(jammerChannels)
        rays_jammers{i} = jammerChannels{i}.VisualVar.Rays;
        if numel(rays_jammers{i}) > 0
            plot(rays_jammers{i}, Type="pathloss", ColorLimits=[40 100]);
        end
    end
end
% Initialize data tables
RSSI_headers = arrayfun(@(x) sprintf('RSSI_Locator_%d', x), 1:numLocators, 'UniformOutput', false);
data = table('Size', [numNodes, 4 + numLocators], ...
             'VariableTypes', [repmat("double", 1, 4 + numLocators)], ...
             'VariableNames', [{'NodeIndex', 'NodeX', 'NodeY', 'NodeZ'}, RSSI_headers]);

data_jammed = table('Size', [numNodes, 4 + numLocators], ...
             'VariableTypes', [repmat("double", 1, 4 + numLocators)], ...
             'VariableNames', [{'NodeIndex', 'NodeX', 'NodeY', 'NodeZ'}, RSSI_headers]);

% Initialize waveforms storage
waveforms = cell(numNodes, numLocators);
waveforms_jammed = cell(numNodes, numLocators);

% Main processing loop
for a = 1:numNodes
    NodePos = Posnodes(:, a)';
    data.NodeIndex(a) = a;
    data.NodeX(a) = NodePos(1);
    data.NodeY(a) = NodePos(2);
    data.NodeZ(a) = NodePos(3);
    data_jammed.NodeIndex(a) = a;
    data_jammed.NodeX(a) = NodePos(1);
    data_jammed.NodeY(a) = NodePos(2);
    data_jammed.NodeZ(a) = NodePos(3);
    RSSI_values = zeros(1, numLocators);  
    RSSI_values_jammed = zeros(1, numLocators);
    
    for i = 1:numLocators
        % Generate transmit bits and waveform
        txBits = randi([0 1], dataLength * bitsPerByte, 1, "int8"); 
        txWaveform = bleWaveformGenerator(txBits, Mode=simMode, ...
                SamplesPerSymbol=sps, ...
                ChannelIndex=rxPhys{a}.ChannelIndex, ...
                AccessAddress=rxPhys{a}.AccessAddress);
        txWaveform = txWaveform * sqrt(P_signal);
        
        % Re-initialize impairments for this iteration
        initImp = initImpTemplate;
        initImp.pfo.FrequencyOffset = randsrc(1,1,-50e3:10:50e3);
        initImp.pfo.PhaseOffset = randsrc(1,1,-10:5:10);
        initoff = 0.15 * sps;
        stepsize = 20 * 1e-6;
        initImp.vdelay = initoff + stepsize * (0:length(txWaveform)-1)';
        initImp.dc = 20;
        txImpairedWfm = helperBLEImpairmentsAddition(txWaveform, initImp);
        
        % *** Added: Reset the fading channel for the signal ***
        reset(rxPhys{a}.BLEChannels{i}.fadingChan);
        
        % Pass through fading channel
        chanDelay = info(rxPhys{a}.BLEChannels{i}.fadingChan).ChannelFilterDelay;
        txChanWfm = rxPhys{a}.BLEChannels{i}.fadingChan([txImpairedWfm; zeros(chanDelay,1)]);
        txChanWfm = txChanWfm(chanDelay+1:end, 1);
        
        % Add AWGN
        rxWaveform = awgn(txChanWfm, snr_dB, "measured");
        rxWaveform = rxWaveform - mean(rxWaveform);
        waveforms{a, i} = rxWaveform;
        RSSI = calculateRSSI(rxWaveform);
        RSSI_values(i) = RSSI;
        
        % Process with jammer
        jammerBits = randi([0 1], dataLength * bitsPerByte, 1, "int8");
        jammerWaveform = bleWaveformGenerator(jammerBits, Mode=simMode, SamplesPerSymbol=sps);
        jammerWaveform = jammerWaveform * sqrt(jammerPower_W);
        % *** Added: Reset the fading channel for the jammer ***
        reset(jammerChannels{a}.fadingChan);
        
        chanDelayJammer = info(jammerChannels{a}.fadingChan).ChannelFilterDelay;
        jammerChanWaveform = jammerChannels{a}.fadingChan([jammerWaveform; zeros(chanDelayJammer, 1)]);
        jammerChanWaveform = jammerChanWaveform(chanDelayJammer + 1:end, 1);
        
        % Ensure waveforms are the same length
        minLength = min(length(txChanWfm), length(jammerChanWaveform));
        txChanWfm = txChanWfm(1:minLength);
        jammerChanWaveform = jammerChanWaveform(1:minLength);
        
        % Combine waveforms and add AWGN
        combinedWaveform = txChanWfm + jammerChanWaveform;
        P_noise = N0 * B + mean(abs(jammerChanWaveform).^2);
        P_signal_rx = mean(abs(txChanWfm).^2);
        SNR = P_signal_rx / P_noise;
        SNR_dB_new = 10 * log10(SNR);
        
        % Ensure SNR_dB_new is within a reasonable range
        SNR_dB_new = max(SNR_dB_new, -20); % Limit minimum SNR to -20 dB
        rxWaveformJammed = awgn(combinedWaveform, SNR_dB_new, 'measured');
        waveforms_jammed{a, i} = rxWaveformJammed;
        RSSI_jammed = calculateRSSI(rxWaveformJammed);
        RSSI_values_jammed(i) = RSSI_jammed;
    end
    data{a, 5:end} = RSSI_values;
    data_jammed{a, 5:end} = RSSI_values_jammed;
end

load('all_random_forest_models.mat');

% Prepare the predictors (RSSI values) for and _jammed
predictors_ = data{:, 5:end};
predictors_jammed = data_jammed{:, 5:end};

% Predict positions using the trained models 
predictedX_ = predict(models{1}, predictors_);
predictedY_ = predict(models{2}, predictors_);
predictedZ_ = predict(models{3}, predictors_);

% Predict positions using the trained models _jammed
predictedX_jammed = predict(models{1}, predictors_jammed);
predictedY_jammed = predict(models{2}, predictors_jammed);
predictedZ_jammed = predict(models{3}, predictors_jammed);

% Extract true positions
trueX = data.NodeX;
trueY = data.NodeY;
trueZ = data.NodeZ;

% Compute errors 
errorsX_ = trueX - predictedX_;
errorsY_ = trueY - predictedY_;
errorsZ_ = trueZ - predictedZ_;

% Compute errors _jammed
errorsX_jammed = trueX - predictedX_jammed;
errorsY_jammed = trueY - predictedY_jammed;
errorsZ_jammed = trueZ - predictedZ_jammed;

% Compute mean absolute errors 
MAE_X_ = mean(abs(errorsX_));
MAE_Y_ = mean(abs(errorsY_));
MAE_Z_ = mean(abs(errorsZ_));

% Compute mean absolute errors _jammed
MAE_X_jammed = mean(abs(errorsX_jammed));
MAE_Y_jammed = mean(abs(errorsY_jammed));
MAE_Z_jammed = mean(abs(errorsZ_jammed));

% Compute overall position errors
positionErrors_ = sqrt(errorsX_.^2 + errorsY_.^2 + errorsZ_.^2);
MAE_position_ = mean(positionErrors_);

positionErrors_jammed = sqrt(errorsX_jammed.^2 + errorsY_jammed.^2 + errorsZ_jammed.^2);
MAE_position_jammed = mean(positionErrors_jammed);

% Display results
fprintf(' Before Jamming Attack:\n');
fprintf('Mean Absolute Error in X: %.4f meters\n', MAE_X_);
fprintf('Mean Absolute Error in Y: %.4f meters\n', MAE_Y_);
fprintf('Mean Absolute Error in Z: %.4f meters\n', MAE_Z_);
fprintf('Mean Absolute Position Error: %.4f meters\n', MAE_position_);

fprintf('\nAfter Jamming Attack:\n');
fprintf('Mean Absolute Error in X: %.4f meters\n', MAE_X_jammed);
fprintf('Mean Absolute Error in Y: %.4f meters\n', MAE_Y_jammed);
fprintf('Mean Absolute Error in Z: %.4f meters\n', MAE_Z_jammed);
fprintf('Mean Absolute Position Error: %.4f meters\n', MAE_position_jammed);

%% Visualization
% Visualize the true positions and predicted positions  and _jammed

% Visualize 
figure;
subplot(1, 2, 1);
title('Before Jamming Attack');
metric_ = visualizePositioningResults('bigeli.stl', [trueX, trueY, trueZ], [predictedX_, predictedY_, predictedZ_]);

% Visualize _jammed
subplot(1, 2, 2);
title('After Jamming Attack');
metric_jammed = visualizePositioningResults('bigeli.stl', [trueX, trueY, trueZ], [predictedX_jammed, predictedY_jammed, predictedZ_jammed]);

% Plot CDF of position errors
figure;
hold on;
ecdf(positionErrors_);
ecdf(positionErrors_jammed);
legend('Normal', 'Jammed');
xlabel('Position Error (m)');
ylabel('Cumulative Probability');
title('CDF of Position Errors');
grid on;
hold off;

% Initialize spectrum analyzer object
spectrumAnalyzerObj = spectrumAnalyzer(...
    'Title', 'Spectrum Before and After Jamming', ...
    'SpectrumType', 'Power density', ...
    'PlotAsTwoSidedSpectrum', true, ...  % Set to true for complex data
    'SampleRate', sampleRate, ...
    'ShowLegend', true);

% Loop through each node and locator to visualize waveforms
for a = 1:numNodes
    for i = 1:numLocators
        % Retrieve the waveforms
        waveform_ = waveforms{a, i};  % *** Modified: Changed waveforms_ to waveforms ***
        waveform_jammed = waveforms_jammed{a, i};
        
        % Ensure waveforms are of the same length
        minLength = min(length(waveform_), length(waveform_jammed));
        waveform_ = waveform_(1:minLength);
        waveform_jammed = waveform_jammed(1:minLength);
        
        % Visualize
        spectrumAnalyzerObj([waveform_, waveform_jammed]);
        pause(0.01);
    end
end

% Function to visualize positioning results
function meanErr = visualizePositioningResults(mapFileName, truthPos, predictedPos)
    % Visualize indoor location estimation predictions
    % truthPos: Nx3 matrix of true positions
    % predictedPos: Nx3 matrix of predicted positions

    % Triangulate the map
    if ~isa(mapFileName, 'triangulation')
        tri = stlread(mapFileName);
    else
        tri = map;
    end
    trisurf(tri, ...
        'FaceAlpha', 0.3, ...
        'FaceColor', [.5, .5, .5], ...
        'EdgeColor', 'none');
    hold on; axis equal; grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    view([84.75 56.38]);

    % Plot edges
    fe = featureEdges(tri, pi/20);
    numEdges = size(fe, 1);
    pts = tri.Points;
    aPts = pts(fe(:, 1), :);
    bPts = pts(fe(:, 2), :);
    fePts = cat(1, reshape(aPts, 1, numEdges, 3), ...
        reshape(bPts, 1, numEdges, 3), nan(1, numEdges, 3));
    fePts = reshape(fePts, [], 3);
    plot3(fePts(:, 1), fePts(:, 2), fePts(:, 3), 'k', 'LineWidth', .5);

    % Compute the distance error
    mErr = sqrt(sum((truthPos - predictedPos).^2, 2));

    % Set the color bar properties
    minErr = floor(min(mErr)) * 10;
    maxErr = ceil(max(mErr)) * 10;
    numColors = (maxErr - minErr) / 5;
    cm = colormap(jet(numColors));

    % Plot the true receiver locations - colored by magnitude of error
    for i = 1:size(truthPos, 1)
        cmIdx = find((mErr(i) * 10 - (minErr:5:maxErr)) < 0, 1) - 1;
        if isempty(cmIdx) || cmIdx < 1
            cmIdx = 1;
        elseif cmIdx > numColors
            cmIdx = numColors;
        end
        scatter3(truthPos(i, 1), truthPos(i, 2), truthPos(i, 3), 50, cm(cmIdx, :), 'filled');
    end

    % Plot predicted positions
    scatter3(predictedPos(:, 1), predictedPos(:, 2), predictedPos(:, 3), 50, 'kx');

    % Create colorbar
    cb = colorbar;
    cb.Label.String = 'Distance error (m)';
    cbLim = cb.Limits;
    cb.Ticks = linspace(cbLim(1), cbLim(2), numColors + 1);
    cb.TickLabels = arrayfun(@(x) sprintf('%.1f', x / 10), linspace(minErr, maxErr, numColors + 1), 'UniformOutput', false);
    title({'{\bf{Position Prediction}}'; '\fontsize{10}Actual Positions Colored by Distance Error'}, 'FontWeight', 'Normal');

    % Compute mean error
    meanErr = mean(mErr);

    hold off;
end

% Target position from the image
targetPosition = [-5.5119, -5.8000, -0.0791]; % Replace with the values from the image

% Find the node index by matching position
[~, nodeIndex] = min(vecnorm(data{:, 2:4} - targetPosition, 2, 2));

% Extract and display RSSI values
rssi_ = data{nodeIndex, 5:end};
rssi_jammed = data_jammed{nodeIndex, 5:end};

fprintf('RSSI Values for Node at Position (%.4f, %.4f, %.4f):\n', targetPosition);
disp(rssi_);
fprintf('RSSI Values After Jamming for Node at Position (%.4f, %.4f, %.4f):\n', targetPosition);
disp(rssi_jammed);
