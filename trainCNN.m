% Load the data
rng(1);
data = readtable('RSSI_data.csv');

% Separate the response (target) and predictors
response = data{:, 1:3};          % First three columns as the target variables
predictors = data{:, 4:end};      % All columns from the fourth onward as predictors

% Split the data into training and test sets
cv = cvpartition(size(data,1), 'HoldOut', 0.2);  % 80% training, 20% testing
trainIdx = cv.training;
testIdx = cv.test;

% Create training and test sets
trainPredictors = predictors(trainIdx, :);
testPredictors = predictors(testIdx, :);
trainResponse = response(trainIdx, :);
testResponse = response(testIdx, :);

% Get the number of features
numFeatures = size(trainPredictors,2);
numResponses = size(trainResponse,2);

% Normalize the predictors
predictorMeans = mean(trainPredictors, 1);
predictorStds = std(trainPredictors, 0, 1);
trainPredictorsNorm = (trainPredictors - predictorMeans) ./ predictorStds;
testPredictorsNorm = (testPredictors - predictorMeans) ./ predictorStds;

% Normalize the responses
responseMeans = mean(trainResponse, 1);
responseStds = std(trainResponse, 0, 1);
trainResponseNorm = (trainResponse - responseMeans) ./ responseStds;
testResponseNorm = (testResponse - responseMeans) ./ responseStds;

% Reshape predictors for CNN input
trainPredictorsReshaped = reshape(trainPredictorsNorm', [numFeatures,1,1,size(trainPredictorsNorm,1)]);
testPredictorsReshaped = reshape(testPredictorsNorm', [numFeatures,1,1,size(testPredictorsNorm,1)]);

% Define the CNN layers
layers = [
    imageInputLayer([numFeatures 1 1], 'Name', 'input')
    convolution2dLayer([5 1], 16, 'Padding', 'same', 'Name', 'conv1')
    batchNormalizationLayer('Name', 'bn1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(64, 'Name', 'fc1')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(numResponses, 'Name', 'fc2')
    regressionLayer('Name', 'output')];

% Specify the training options
options = trainingOptions('adam', ...
    'MiniBatchSize', 32, ...
    'MaxEpochs', 30, ...
    'InitialLearnRate',1e-3, ...
    'Plots','training-progress', ...
    'Verbose',false);

% Train the CNN model
net = trainNetwork(trainPredictorsReshaped, trainResponseNorm, layers, options);

% Predict on the test data
predictedResponseNorm = predict(net, testPredictorsReshaped);

% Unnormalize the predicted responses
predictedResponse = predictedResponseNorm .* responseStds + responseMeans;

% Calculate mean absolute error (MAE) on the test set
mae_test = mean(abs(predictedResponse - testResponse));

% Display MAE for each target variable
for i = 1:numResponses
    fprintf('Mean Absolute Error for Target %d (Test Set): %.4f\n', i, mae_test(i));
end

% Save the model and normalization parameters
save('cnn_model.mat', 'net', 'predictorMeans', 'predictorStds', 'responseMeans', 'responseStds', '-v7.3');
