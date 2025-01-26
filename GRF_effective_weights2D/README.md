<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: GRF_effective_weights2D

Published in: GRF

Description: Estimation of effective weights alpha_i(x_1, x_2) for the infeasable observations theta_tilde. The effective weights are computed by a regression forest on a grid of observations x_ij=(-0.5 + i/n1, 0 + 0.02 * j) for i=1,...,n1 and j=0,...,50, for n1=50, 100, 200 and target variable Y_i=theta(x_ij) + eps_ij, for a given theta function, here triangle function theta(x_ij) = max(0, 1 - |x_ij,1|/0.2), with Gaussian noise eps_ij with mean zero and standard deviation of 0 and 0.1.

Keywords: RF, GRF, infeasable function, estimation, effective weights, bandwidth, contour

Author: Marius Sterling

See also: 

Submitted: 20.01.2021

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_effective_weights2D/RF_theta_triangle___effective_weights___sigma000___n050___contour___x_000_050.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_effective_weights2D/RF_theta_triangle___effective_weights___sigma000___n100___contour___x_000_050.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_effective_weights2D/RF_theta_triangle___effective_weights___sigma000___n200___contour___x_000_050.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_effective_weights2D/RF_theta_triangle___effective_weights___sigma010___n050___contour___x_000_050.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_effective_weights2D/RF_theta_triangle___effective_weights___sigma010___n100___contour___x_000_050.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/GRF/master/GRF_effective_weights2D/RF_theta_triangle___effective_weights___sigma010___n200___contour___x_000_050.png" alt="Image" />
</div>

