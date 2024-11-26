<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Published in: 'GRF'

Description: 'Estimation of effective weights alpha_i(x_1) for the infeasable observations theta_tilde. The effective weights are computed by a quantile regression forest on a grid of observation for given n and target variable Y_i=theta(x_ij) + eps_ij, for a given theta function, here triangle function theta(x_ij) = max(0, 1 - |x_ij,1|/0.2), with Gaussian noise eps_ij with mean zero and standard deviation of 0 and 0.1.' The theta is then estimated and the estimated function is tested on a new grid to calculate coverage of true theta by the point wise confidence intervals created using multiplier bootstrap. At last, simultaneous confidence intervals are calculated and plotted

Keywords: 'RF, GRF, infeasable function, estimation, effective weights, confidence bands, confidence intervals '

See also: ''

Author: 'Kainat Khowaja'

Submitted: '10.05.2021'

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_confidence_bands_1D/CI_bands_n5000_tau005_sig001_grids100_.png" alt="Histogram" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_confidence_bands_1D/CI_bands_n5000_tau009_sig001_grids100_.png" alt="Histogram" />
</div>

