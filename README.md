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

- time range can be calculated autoatically or keep fixed (see options)
- 2D histogram added: x-axis: time y-axis: pulse height (similar to time slicing)