% Load the data
rng(1);
data = readtable('RSSI_data.csv');
% Separate the response (target) and predictors
response = data{:, 1:3};          % First three columns as the target variables
predictors = data(:, 4:end);      % All columns from the fourth onward as predictors

% Split the data into training and test sets
cv = cvpartition(height(data), 'HoldOut', 0.2);  % 80% training, 20% testing
trainIdx = cv.training;
testIdx = cv.test;

% Create training and test sets
trainPredictors = predictors(trainIdx, :);
testPredictors = predictors(testIdx, :);
trainResponse = response(trainIdx, :);
testResponse = response(testIdx, :);

% Initialize parameters
numTrees = 200;
models = cell(1, 3);  % Cell array to store each model

% Train a random forest model for each target variable on training data
for i = 1:3
    models{i} = TreeBagger(numTrees, trainPredictors, trainResponse(:, i), 'Method', 'regression');
end

% Initialize an array to store mean absolute errors for test set
mae_test = zeros(1, 3);

% Calculate the mean absolute error for each target variable on test set
for i = 1:3
    % Predict on the test data for each target variable
    predictedLabels = predict(models{i}, testPredictors);
    
    % Calculate mean absolute error (MAE) on the test set
    mae_test(i) = mean(abs(predictedLabels - testResponse(:, i)));
    fprintf('Mean Absolute Error for Target %d (Test Set): %.4f\n', i, mae_test(i));
end

% Save all models to a single file
save('all_random_forest_models.mat', 'models', '-v7.3');
