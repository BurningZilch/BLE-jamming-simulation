# BLE-simulation
This in an assessment project in my university that aims to simulate how Bluetooth Low Energy (BLE) interference impacts the performance of a machine learning model. Due to the time pressure during the development, there might be error in the code, please use with caution.

## **How to Run the Project**

1. **Step 1: Simulate BLE Data**
   - Run `generate_dataset.m` to create the required training dataset.
   - Modify `createScenario.m` and `createJammer.m` to simulate specific BLE environments.

2. **Step 2: Train Models**
   - Run `trainCNN.m` to train the CNN model.
   - Run `trainRF.m` to train the RF model.

3. **Step 3: Test Models**
   - Use `testCNN.m` to evaluate the CNN model.
   - Use `testRandomForests.m` to evaluate the RF model.

