# BLE-Simulation

This project is part of a university assessment that simulates how **Bluetooth Low Energy (BLE)** jamming interference impacts the performance of a machine learning model. Due to time constraints during development, there may be some errors in the code. **Please use it with caution.**

---

## Real-World Scenario  
In practice, BLE devices typically send **10–100 signals per second**, and devices calculate an **average value** from these signals for location purposes.  
⚠️ **Note:** This behavior is **not implemented** in the current code.

---
## Website
This project is accompanied by a blog that documents the development process. You can read it here: [https://burningzilch.github.io/2024-11-25-BLEMutipathFading/](https://burningzilch.github.io/2024-11-25-BLEMutipathFading/)

---
## Issues and Questions  
If you encounter any issues or have questions about this project, feel free to reach out.


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



