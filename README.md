# FloodHazardModelCalibrationUsingMultiresolutionModelOutput
https://arxiv.org/abs/2203.00840

Base folder:  .../FloodingModelCalibrationProject/multires/spatialchangPC2014/10m50m/code/nCh_RWE/

To hold out extremes: 
- Folder: Base folder
- File to hold out extremes in 50e200c combo: holdOutExtremes50e200c.R
- File to hold out extremes in 100e400c combo: holdOutExtremes100e400c.R

To use 100 expensive model runs and 400 cheap model runs, go to folder: 100e400c
To use 50 expensive model runs and 200 cheap model runs, go to folder: 50e200c

After selecting the desired combination of model runs:

- Folder: base folder + 100e400c/CompareCalibration

To compare RMSEs among emulators with different combinations of model runs: 
- Folder : Base folder
- File: compareRMSEs.R

To compare calibration times between approaches with different combinations of model runs: 
- Folder : Base folder
- File: compareTimes.R
