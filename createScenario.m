function [locators, nodes, antPosNode] = createScenario(Posnodes,Poslocators)
%createScenario Generates the tx locator and rx node objects
%based on separation between the nodes

fc = 2.4e9; % Set the carrier frequency (Hz) to one of broadcasting channels
lambda = physconst("lightspeed")/fc;

txArray = arrayConfig("Size", [1 1], "ElementSpacing", 2*lambda);
rxArray = arrayConfig("Size", [1 1], "ElementSpacing", lambda);

% Update `antPosLocator` and `antPosNode` with these arrays
antPosLocator = Poslocators;
antPosNode = Posnodes;

% Create multiple locator and node sites with a single constructor call
locators = txsite("cartesian", ...
    AntennaPosition=antPosLocator, ...
    Antenna=txArray, ...
    TransmitterFrequency=fc);

nodes = rxsite("cartesian", ...
    AntennaPosition=antPosNode, ...
    Antenna=rxArray, ...
    AntennaAngle=[0;90]);

end