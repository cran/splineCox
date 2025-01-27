# splineCox 0.0.1

- Initial release on CRAN.
- Implements a two-stage estimation approach for Cox regression.
- Uses five-parameter M-spline functions for baseline hazard modeling.
- Provides model selection based on log-likelihood criteria.

* Initial CRAN submission.

# splineCox 0.0.2

- Fixed bugs in the `splineCox.reg2` function when handling custom numeric vectors.
- Added examples for custom vector normalization in the vignette.
- Improved documentation and error messages.


# splineCox 0.0.3

## New Features
- Added plotting functionality to `splineCox.reg1` and `splineCox.reg2` for visualizing the estimated baseline hazard function with confidence intervals.
- Updated the vignette to include descriptions and examples of the plotting feature.

## Documentation Improvements
- Improved documentation for `splineCox.reg1` and `splineCox.reg2` to clearly describe the new plotting feature.
- Updated the vignette to demonstrate how to use the plotting functionality.
