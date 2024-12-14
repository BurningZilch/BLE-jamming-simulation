function rssi = calculateRSSI(waveform)
    dBdBmConvFactor = 30;
    btWaveformPower = var(waveform);
    rssi = 10 * log10(btWaveformPower) + dBdBmConvFactor;
end