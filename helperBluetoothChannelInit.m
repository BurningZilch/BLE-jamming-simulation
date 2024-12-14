function init = helperBluetoothChannelInit(sampleRate,channelModel,txPosition,rxPosition)
% helperBluetoothChannelInit initializes multipath fading channel model for
% Bluetooth waveforms in indoor scenarios.
%
%   INIT = helperBluetoothChannelInit(SAMPLERATE,CHANNELMODEL) outputs a
%   structure, INIT, which contains the front-end multi-path fading channel
%   parameters for the given SAMPLERATE and CHANNELMODEL.
%
%   SAMPLERATE is a scalar quantity specifying the sampling frequency of
%   the Bluetooth waveform.
%
%   CHANNELMODEL is a character vector or string specifying the type of
%   channel model. It must be one of the following: "RAYLEIGH CHANNEL",
%   "RICIAN CHANNEL", and "RAYTRACING CHANNEL".
%
%   References: [1] International Telecommunication Union, ITU-M
%   Recommendation M.1225 GUIDELINES FOR EVALUATION OF RADIO TRANSMISSION
%   TECHNOLOGIES FOR IMT-2000, Annex 2, Section 1.2.2 Channel impulse
%   response model, Table 3 "Indoor office test environment
%   tapped-delay-line parameters".
%
%   See also bleWaveformGenerator, bluetoothWaveformGenerator.

%   Copyright 2022-2024 The MathWorks, Inc.

% Create multipath fading channel object
if any(channelModel == ["Rician Channel","Rayleigh Channel"])
    % Define channel delay profile for channel model type "Rician Channel"
    % and "Rayleigh Channel"

    % Select the channel delay profile as "Channel A" or "Channel B", as
    % specified in [1].
    channelDelayProfile = "Channel A";

    % Select the taps for the channel model, as a scalar or a row
    % vector. When you set tapNumber to a scalar, the channel is frequency
    % flat. When you set tapNumber to a vector, the channel is frequency
    % selective.
    tapNumber = [3 6];

    % Specify relative delays and average power for channel A and channel
    % B.
    relativeDelaysChannelA = [0 50 110 170 290 310]*1e-9;   % In seconds
    avgPathGainsChannelA = [0 -3 -10 -18 -26 -32];          % in dB
    relativeDelaysChannelB = [0 100 200 300 500 700]*1e-9;  % In seconds
    avgPathGainsChannelB = [0 -3.6 -7.2 -10.8 -18 -25.2];   % in dB

    % Based on the type of channel and number of taps specified, calculate
    % path delays and average path gains for the channel model.
    if channelDelayProfile == "Channel A"
        pathDelays = relativeDelaysChannelA(tapNumber);
        avgPathGains = avgPathGainsChannelA(tapNumber);
    elseif channelDelayProfile == "Channel B"
        pathDelays = relativeDelaysChannelB(tapNumber);
        avgPathGains = avgPathGainsChannelB(tapNumber);
    end

    % Create fading channel object
    switch channelModel
        case "Rayleigh Channel"
            % Initialize Rayleigh channel
            init.fadingChan = comm.RayleighChannel(SampleRate=sampleRate, ...
                RandomStream="Global stream", ...
                AveragePathGains=avgPathGains, ...
                PathDelays=pathDelays);

        case "Rician Channel"
            % Initialize Rician channel
            init.fadingChan = comm.RicianChannel(SampleRate=sampleRate, ...
                RandomStream="Global stream", ...
                AveragePathGains=avgPathGains, ...
                PathDelays=pathDelays);
    end
elseif channelModel=="Raytracing Channel"
    % Define conference room environment for raytracing
    mapFileName = "bigeli.stl";

    % Define transmitter site and receiver site
    tx = txsite("cartesian", ...
        AntennaPosition=txPosition, ...
        TransmitterFrequency=2.426e9);
    rx = rxsite("cartesian", ...
        AntennaPosition=rxPosition);

    % Define propagation model
    pm = propagationModel("raytracing", ...
        CoordinateSystem="cartesian", ...
        MaxNumReflections=2, ...
        SurfaceMaterial="wood");

    % Plot rays in the conference room
    rays = raytrace(tx,rx,pm,Map=mapFileName,Type="pathloss");
    
    if numel(rays{1}) > 0
      
        rays = rays{1,1}; 

    % Initialize Raytracing channel
        init.fadingChan = comm.RayTracingChannel(rays,tx,rx);
        init.fadingChan.SampleRate = sampleRate;
        init.fadingChan.ReceiverVirtualVelocity = randsrc(3,1,0:0.001:0.01);
    else %TODO something is wroing here
        init.fadingChan = comm.RayleighChannel(SampleRate=sampleRate, ...
                                           PathDelays=0, ...
                                           AveragePathGains=0);

    end
    init.VisualVar = struct("MapFileName",mapFileName,"TxSite",tx,"RxSite",rx,"Rays",rays);
   
end
end