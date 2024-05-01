# Risk-Averse Optimization-based Outlier Accommodation for RTK INS Fusion
UrbanRTK-INS-OutlierOpt is an open-source framework that integrates Real-Time Kinematic (RTK) 
Global Navigation Satellite Systems (GNSS) with Inertial Navigation Systems (INS), 
utilizing a Risk-Averse Performance-Specified (RAPS) Optimization approach. 
This repository is designed to provide robust outlier accommodation in urban environments,
where GNSS signals are often compromised due to obstacles like buildings and bridges.

# RAPS-RTK-INS Framework
![BlockDiagram_Updated](https://github.com/Azurehappen/UrbanRTK-INS-OutlierOpt/assets/45580484/3fbf5612-53aa-4845-ba98-f3f8237f764f)

# Experimental Results
The open-source [TEX-CUP](https://radionavlab.ae.utexas.edu/texcup-desc/) dataset (2019May09) is used.  The experimental route traversing areas within the west campus of The University of Texas at Austin and downtown Austin, contains viaducts, high-rise buildings, and dense foliage. Results are estimated through forward (real-time) processing.

The RTK-GNSS/INS integration utilizes single-frequency measurements from a Septentrio receiver GPS L1, GLONASS L1, GALILEO E1, and Beidou B1.  The inertial measurements are provided by the LORD MicroStrain 3DM GX5-25 IMU with a sampling rate of 100 Hz.

## Full trajectory from RTK-INS-RAPS
![image](https://github.com/Azurehappen/UrbanRTK-INS-OutlierOpt/assets/45580484/2dfd3020-187d-4b78-928b-36fad6dad2c1)

## Detailed Comparison

Left panels: near the Dell Medical School buildings.

Right panels: near Sailboat Building (multiple skyscrapers surround the site)

### Results from the traditional Threshold Decision (TD) method (TD-RTK-INS)
![td_area](https://github.com/Azurehappen/UrbanRTK-INS-OutlierOpt/assets/45580484/7cb2d67f-8353-44cf-a4d1-294541f425fc)

### Results from RAPS-RTK-INS
![raps_area](https://github.com/Azurehappen/UrbanRTK-INS-OutlierOpt/assets/45580484/40a5d462-f4a1-412f-9e7b-7319615f263b)
