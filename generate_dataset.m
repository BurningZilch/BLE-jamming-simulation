rng(1);
%init 
k = 1.38e-23;  % Boltzmann constant in J/K
T = 290;       % Temperature in Kelvin
%for calculate thermal noise

% Parameters
P_signal = 0.001;    % 1 mW
B = 2e6;           % Bandwidth in Hz 
R_b = 125e3;         % Bit rate in bits per second 
P_dBm = 10 * log10(P_signal / 0.001);
% Calculate Noise Density N0
N0 = k * T;           % Noise density in W/Hz
N0_dBmHz = 10 * log10(N0) + 30;  % Noise density in dBm/Hz
snr_dB = 20;
% Calculate Noise Power
P_noise = N0 * B;

fc = 2.426e9;               % channel 38 
c = 3e8;                  % Speed of light (m/s)
lambda = c / fc;          % Wavelength (m)

sps = 2;
sampleRate = sps*125e3; 
dataLength = 128;         % Includes header, payload, and CRC
phyMode = "BR";
simMode = "LE1M";
%packetType = "DM1";
channelModel = "Raytracing Channel";
enableEqualizer = false;
bitsPerByte = 8;
%--------------------------------------------------------------
%node generate
nodeDistance = 2; 
xRange = [-5.6, 5.6]; 
yRange = [-5.6, 5.6]; 
zRange = [0, 4]; 


% Generate positions
Posnodes = generatePosnodes(nodeDistance, xRange, yRange, zRange);
disp(Posnodes)

AdvertisingAddress = [0 1 1 0 1 0 1 1 0 1 1 1 1 1 0 1 1 0 0 1 0 0 0 1 0 1 1 1 0 0 0 1]';

Poslocators = [
    2.8, 3.0, -1.8, 5.0, -2.2, 1.4, 4.0, -2.8, 3.6, -1.2;
   -2.0, 4.2,  6.0, -1.6, -4.4, 3.0, 0.0, -3.4, 4.6,  2.4;
    4.0, 4.0, 2.4, 1.2, 3.6, -0.8, 6.6, -4.0, -1.0, 5.4];

[locators,nodes,antPosNode] = createScenario(Posnodes,Poslocators);
initImp = helperBLEImpairmentsInit("LE Coded PHY",sps);

spectrumAnalyzerBasic = spectrumAnalyzer(...
    'Title', 'Spectrum Before and After Jamming', ...
    'SpectrumType', 'Power density', ...
    'PlotAsTwoSidedSpectrum', true, ...  % Set to true for complex data
    'SampleRate', sampleRate, ...
    'ShowLegend', true);

%make PHY object---------------------------------------
rxPhys = cell(1,1);
txPhys = cell(1,1);
for i = 1:size(Poslocators, 2)
    txPhy = helperBluetoothPHY('TxPower', P_dBm, 'Mode', phyMode);
    txPhy.NodePosition = locators(i).AntennaPosition;
    txPhys{i} = txPhy;
    
end
for i = 1:size(Posnodes, 2)
    rxPhy = helperBluetoothPHY('Mode', phyMode);
    rxPhy.NodePosition = nodes(i).AntennaPosition;
    rxPhys{i} = rxPhy;
end
% Channel
for i = 1:size(Posnodes, 2) 
    for a = 1:size(Poslocators, 2) 
        BLEChannel = helperBluetoothChannelInit(sampleRate,channelModel,locators(a).AntennaPosition,nodes(i).AntennaPosition);
        rxPhys{i}.BLEChannels{a} = BLEChannel;
    end
end
%show the room-------------------------------------------
 viewer = siteviewer(SceneModel=BLEChannel.VisualVar.MapFileName);
% 
 show(nodes,Icon="bleRxIcon.png");
 show(locators,Icon="bleTxIcon.png");
 nodenum = 2; % which node to display

 rays = cell(1,1);
 for i = 1:length(rxPhys{nodenum}.BLEChannels)
     rays{i} = rxPhys{nodenum}.BLEChannels{i}.VisualVar.Rays;
     if numel(rays{i}) > 0
     plot(rays{i}, Type="pathloss", ColorLimits=[0 120]);
     end
 end
%------------------------------------------------------------
%generate data for packet, and get wave from data
numNodes = size(Posnodes, 2);
numLocators = size(Poslocators, 2);
RSSI_headers = arrayfun(@(x) sprintf('RSSI_Locator_%d', x), 1:numLocators, 'UniformOutput', false);
data = table('Size', [numNodes, 4 + numLocators], ...
             'VariableTypes', [repmat("double", 1, 4 + numLocators)], ...
             'VariableNames', [{'NodeIndex', 'NodeX', 'NodeY', 'NodeZ'}, RSSI_headers]);
numDataSets = 32; % make 32 datasets for each node 
allData = table('Size', [0, 4 + numLocators], ...
                'VariableTypes', [repmat("double", 1, 4 + numLocators)], ...
                'VariableNames', [{'NodeIndex', 'NodeX', 'NodeY', 'NodeZ'}, RSSI_headers]);

for dSet = 1:numDataSets
for a = 1:size(Posnodes, 2)  % for node a
 NodePos = Posnodes(:, a)';
    data.NodeIndex(a) = a;
    data.NodeX(a) = NodePos(1);
    data.NodeY(a) = NodePos(2);
    data.NodeZ(a) = NodePos(3);
    RSSI_values = zeros(1, numLocators);
for i = 1:size(Poslocators, 2) % with locator i

txBits = randi([0 1],dataLength*bitsPerByte,1,"int8"); 

txWaveform = bleWaveformGenerator(txBits,Mode=simMode, ...
                SamplesPerSymbol=sps, ...
                ChannelIndex=rxPhys{a}.ChannelIndex, ...
                AccessAddress=rxPhys{a}.AccessAddress);
txWaveform = txWaveform .* sqrt(txWaveform);
%generate RF impairment
initImp.pfo.FrequencyOffset = randsrc(1,1,-50e3:10:50e3);                       % Frequency offset in Hz, range is [-150000,+150000]
initImp.pfo.PhaseOffset = randsrc(1,1,-10:5:10);                                % Phase offset in degrees
initoff = 0.15*sps;                                                             % Static timing offset
stepsize = 20*1e-6;                                                             % Timing drift in ppm, max range is +/- 50 ppm
initImp.vdelay = (initoff:stepsize:initoff+stepsize*(length(txWaveform)-1))';   % Variable timing offset
initImp.dc = 20;                                                                % Percentage related to maximum amplitude value
txImpairedWfm = helperBLEImpairmentsAddition(txWaveform,initImp); % tx which hurts by RF

%Pass the impaired waveform through the fading channel. (multipath)
chanDelay = info(rxPhys{a}.BLEChannels{i}.fadingChan).ChannelFilterDelay;

% Pass through the fading channel model
txChanWfm = rxPhys{a}.BLEChannels{i}.fadingChan([txImpairedWfm; zeros(chanDelay,1)]);
txChanWfm = txChanWfm(chanDelay+1:end,1);

%--------------------------------------------------------------------------
% Add AWGN before receive
rxWaveform = awgn(txChanWfm,snr_dB,"measured");
RSSI = calculateRSSI(rxWaveform);
RSSI_values(i) = RSSI;

% show waveform, we likely do not need this at here
%spectrumAnalyzerBasic(rxWaveform,txWaveform);
%pause(0.01);
end
data{a, 5:end} = RSSI_values;
end
allData = [allData; data];
end
allData = allData(:, 2:end);
%--------------------------------------------------------------------------
%dataset
writetable(allData, 'RSSI_data.csv');

% ------------------------------------------------------------
% Print the first node's RSSI with all locators for debug
fprintf('First node RSSI with all locators:\n');
firstNodeRSSI = allData{1, 4:end};  % Extract RSSI columns for the first node
rssiTable = array2table(firstNodeRSSI, ...
    'VariableNames', allData.Properties.VariableNames(4:end));  % Create a table for better readability
disp(rssiTable);
fprintf('Location of the first node:\n');
firstNodeLocation = allData{1, 1:3};  % Extract X, Y, Z coordinates for the first node
fprintf('X: %.2f, Y: %.2f, Z: %.2f\n', firstNodeLocation(1), firstNodeLocation(2), firstNodeLocation(3));
