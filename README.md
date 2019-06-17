# ESR-Analyses
MATLAB scripts to analyze cw-ESR data. ESR-Analyses requires the
[natural constants](https://github.com/OE-FET/Natural-constants) package.

If you publish any data processed with the ESR-Analyses routines, please cite "Schott, S.
et al. Nat. Phys. (2019) doi:10.1038/s41567-019-0538-0" where the methods implemented here
have been first published.

## About
ESR-Analyses is composed of:

1. General utility functions which are useful in an ESR context:

    - Functions for common resonance lineshapes: `lorentzian`, `gaussian`, etc.
    - Utility functions for common conversions: `b2g` (converts magnetic field to
      g-factor), `chi2nspin` (converts susceptibility to number of spins), etc.
    - Functions to simulate ESR spectra: `field_mod_sim`, `ESRLorentzSimulation`, etc.

2. Functions to read and manipulate Bruker Xepr data files:

    - `BrukerRead` to read Xepr data files and return the measurement data as well all
       measurement parameters.
    - Functions to process the data: `normalize_spectrum`, `subtract_background`,
      `baseline_corr`, etc.

3. Functions to analyse cw-ESR data:

    - Low-level functions for specific tasks: `gfactor_determination`, `double_int_num`,
      `spin_counting`, etc.
    - High-level functions: `PowerSatAnalysesLorentzFit`, `PowerSatAnalysesVoigtFit`, etc.

All functions do exactly what you would expect from their name, and most of them are well
documented. Therefore, please refer to the individual doc-strings for more information.
