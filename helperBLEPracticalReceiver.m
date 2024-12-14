function [bits,accessAddress] = helperBLEPracticalReceiver(rxWaveform,rxCfg)
%helperBLEPracticalReceiver Demodulate and decodes the received Bluetooth
%LE waveform
%
%   [BITS,ACCESSADDRESS] = helperBLEPracticalReceiver(RXWAVEFORM,RXCFG)
%   decodes the received Bluetooth LE waveform,RXWAVEFORM.
%
%   BITS is an int8 column vector containing the recovered information bits
%   with maximum length of 2080 bits.
%
%   ACCESSADDRESS is an int8 column vector of length 32 bits containing the
%   access address information.
%
%   RXWAVEFORM is a complex valued time-domain waveform with size Ns-by-1,
%   where Ns represents the number of received samples.
%
%   RXCFG is a structure containing these fields:
%   "Mode"                  Specify the physical layer reception mode as
%                           one of "LE1M","LE2M","LE500K",and "LE125K".
%
%   "ChannelIndex"          Specify the channel index as an integer in
%                           the range [0,39]. For data channels, specify
%                           this value in the range [0,36]. For advertising
%                           channels, specify this value in the range
%                           [37,39].
%
%   "SamplesPerSymbol"      Specify the samples per symbol as a positive
%                           integer.
%
%   "DFPacketType"          Specify the direction finding packet type as
%                           "ConnectionlessCTE", "ConnectionCTE", or
%                           "Disabled".
%
%   "CoarseFreqCompensator" Specify the coarse frequency compensator system
%                           object as comm.CoarseFrequencyCompensator.
%
%   "PreambleDetector"      Specify the preamble detector system object as
%                           comm.PreambleDetector.
%
%   "EqualizerFlag"         Specify the flag to enable or disable the
%                           equalizer

%   Copyright 2018-2024 The MathWorks, Inc.

% DC offset correction
rxDCFree = rxWaveform - mean(rxWaveform);

% Estimate and compensate frequency offset
rxFreqComp = rxCfg.CoarseFreqCompensator(rxDCFree);
release(rxCfg.CoarseFreqCompensator);

% Generate reference signals used for packet detection
sps = rxCfg.SamplesPerSymbol;
preamble = ble.internal.preambleGenerator(rxCfg.Mode,rxCfg.AccessAddress);
if any(rxCfg.Mode==["LE1M","LE2M"])
    refSequence = [preamble; rxCfg.AccessAddress];
else
    trellis = poly2trellis(4,[17 13]);
    fecAA = convenc(rxCfg.AccessAddress,trellis);
    pattern = [1 1 0 0].';
    patternLen = length(pattern);
    repBlock = reshape(repmat(fecAA.',patternLen,1),1,[]);
    repPattern = reshape(repmat(pattern,1,length(fecAA)),1,[]);
    codedAA = ~xor(repBlock,repPattern).';
    refSequence = [preamble; codedAA];
end
refSamples = ble.internal.gmskmod(refSequence,sps);
if isa(rxWaveform,"single")
    refSamples = single(refSamples);
end

% Perform timing synchronization
prbDet = rxCfg.PreambleDetector;
prbDet.Preamble = refSamples;
[~,dtMt] = prbDet(rxFreqComp);
release(prbDet)
prbDet.Threshold = max(dtMt);
prbIdx = prbDet(rxFreqComp);
release(prbDet)
if prbIdx >= length(prbDet.Preamble)
    rcvTrim = rxFreqComp(1+prbIdx-length(prbDet.Preamble):end);
else
    rcvTrim = rxFreqComp;
end
if rem(length(rcvTrim),sps)
    rcvTrim = [rcvTrim;zeros(sps-rem(length(rcvTrim),sps),1)];
end

% Equalize the fading effects
if any(fieldnames(rxCfg)=="EqualizerFlag") && rxCfg.EqualizerFlag
    % Constellation values for the Bluetooth LE GFSK modulation scheme
    constellationValues = uniquetol([real(refSamples) imag(refSamples)],0.01,"ByRows",1);
    constellationValues = constellationValues(:,1) + 1i*constellationValues(:,2);

    % Define variables required for the decision feedback equalizer
    nForwardFeedTaps = 5;
    nFeedbackTaps = 1;
    refTap = 4;
    stepSize = 1e-5;

    % Define the equalizer
    dfeEq = comm.DecisionFeedbackEqualizer(Algorithm="LMS", ...
        Constellation=constellationValues, ...
        NumForwardTaps=nForwardFeedTaps, ...
        NumFeedbackTaps=nFeedbackTaps, ...
        ReferenceTap=refTap, ...
        StepSize=stepSize);

    % Calculate the DFE delay
    equalizerDelay = info(dfeEq).Latency;

    % Equalization using DFE with the least mean square algorithm
    equalizedWaveform = dfeEq([rcvTrim;zeros(equalizerDelay,1)],refSamples);
    equalizedWaveform = equalizedWaveform(equalizerDelay+1:end,1);
    rcvTrim  = equalizedWaveform;
end

% Recover the data bits
fecBlock1Length = length(refSamples)+40*any(rxCfg.Mode==["LE500K","LE125K"])*sps;
cteInfoInd = 24*(rxCfg.DFPacketType=="ConnectionCTE")+40*(rxCfg.DFPacketType=="ConnectionlessCTE");
crcLength = 24*(any(rxCfg.Mode==["LE1M","LE2M"]))+2*(rxCfg.Mode=="LE500K")+8*(rxCfg.Mode=="LE125K");
term2Length = 3*(2*(rxCfg.Mode=="LE500K")+8*(rxCfg.Mode=="LE125K"));
minPacketLength = fecBlock1Length+(cteInfoInd+crcLength+term2Length)*sps;
if length(rcvTrim) > minPacketLength
    if (rxCfg.Mode=="LE500K") && (rem(length(rcvTrim),2*sps) ~= 0)
        padLen = 2*sps - rem(length(rcvTrim),2*sps);
    elseif (rxCfg.Mode=="LE125K") && (rem(length(rcvTrim),8*sps) ~= 0)
        padLen = 8*sps - rem(length(rcvTrim),8*sps);
    else
        padLen = 0;
    end
    rcvTrim = [rcvTrim;zeros(padLen,1)];
    [bits,accessAddress] = bleIdealReceiver(rcvTrim, ...
        Mode=rxCfg.Mode, ...
        ChannelIndex=rxCfg.ChannelIndex, ...
        SamplesPerSymbol=sps, ...
        AccessAddress=rxCfg.AccessAddress, ...
        DFPacketType=rxCfg.DFPacketType);
else
    bits = [];
    accessAddress = [];
end
end