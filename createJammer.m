function [jammers, antPosJammer] = createJammer(Posjammers)
%createJammer Generates the jammer objects
%based on separation between the jammer positions

fc = 2.4e9; % Set the carrier frequency (Hz) for the jammers
lambda = physconst("lightspeed")/fc;

% Configure the jammer antenna arrays
jammerArray = arrayConfig("Size", [1 1], "ElementSpacing", 1.5*lambda);

% Update `antPosJammer` with the jammer positions
antPosJammer = Posjammers;

% Create multiple jammer sites with a single constructor call
jammers = txsite("cartesian", ...
    AntennaPosition=antPosJammer, ...
    Antenna=jammerArray, ...
    TransmitterFrequency=fc);

end
