# DigiplotQt5

Initial (w/o peak fitting) data analysis tool set,
written by Dr. Boeglin, with minor modifications by A.Netepenko. This tool is to analyze digitizer files produced by NI and GaGe digitizers. Functionality includes:

- display the raw data in a thinned version
- find peaks
- histogram peak heights
- time slice time series
- histogram time sliced
- calculate rates
- create rate profiles
- Basic FFT capabilities

This an ongoing development.  

### v1.0.0 new features:

- read npz files produced from filtering data
- fit histogram peak positions

### v2.0.0 new features:

- time range can be calculated automatically or keep fixed (see options)
- 2D histogram added: x-axis: time,  y-axis: pulse height (similar to time slicing)

### v3.0.0 new features:
 
- using a data channel for measuring electronic background a fitted cross correlated background signal can be subtracted. To view the corrected data select "Corrected" in the options menu. The parameters for this corrections are in the new "Bkg Subtraction Parameters" menu.
- A moving average can be calculated to smooth out high frequency noise. The number of points used for the window size is selected in the "Parameters" menu under "Moving Average".
- A time offset can be added to align the digitizer time with the shot time. This is also in the parameter section under "Time offset". This offset is added to the digitizer time, and can be negative.
- viewing limits for the data can be saved and recalled using "Actions->Save Limits" or  "Actions->Choose Limits". This is usefule to save zoomed areas.
- Found peaks can be cleared using "Actions->Clear Peaks"