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


# splineCox 0.0.3

## New Features
- Added plotting functionality to `splineCox.reg1` and `splineCox.reg2` for visualizing the estimated baseline hazard function with confidence intervals.
- Updated the vignette to include descriptions and examples of the plotting feature.

## Documentation Improvements
- Improved documentation for `splineCox.reg1` and `splineCox.reg2` to clearly describe the new plotting feature.
- Updated the vignette to demonstrate how to use the plotting functionality.


# splineCox 0.0.4

## New Features
- Added a new function `spline.copula` implementing a B-spline copula model based on the five-parameter M-spline basis functions.
- The copula function offers both copula density and distribution function evaluations, with a variety of built-in coefficient matrices (presets) for modeling a wide range of dependence structures.

## Documentation and Publication
- The main methodology of this package, including the spline-based Cox regression has been published in a peer-reviewed journal. 
- Updated documentation to reference the published article.


# splineCox 0.0.5

- Added analytical computation of Kendall’s tau and Spearman’s rho in `spline.copula()`
- Updated documentation and vignette accordingly


# splineCox 0.0.6

- Added spline.copula.simu() to report Kendall's tau and Spearman's rho(empirical and theoretical).

# splineCox 0.0.7
* Fixed missing importFrom(stats, integrate)
* Updated spline.copula.simu with improved functionality