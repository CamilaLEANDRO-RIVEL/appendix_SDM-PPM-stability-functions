# appendix_SDM-PPM-stability-functions
This repository contains the R code to (i) run a species distribution model in the point process modelling framework and to (ii) explore model stability (or the capacity to predict correctly all presence data) within the PPM perspective fitted with a lasso penalty and observer bias corrections. Different R functions which were written to establish intensity and coefficient measures. These functions are the stability assessment toolbox, referred to as “diagnostic tools”. Code, simulated data and a detailed tutorial illustrating use of this code are provided. The model is described in a paper submitted to Ecological Modelling.

Here, a short description of the supplied R functions to explore model stability, which are contained in the DiagnosticFunctions.R file supplied in the supplementary material. Full details of these functions and a demonstration of their usage appears in the RMarkdown tutorial in the supplementary material. 
<table style="width: 826.933px;">
<tbody>
<tr>
<td style="width: 137px;">
<p><strong>Function</strong></p>
</td>
<td style="width: 124px;">
<p><strong>Characteristic</strong></p>
</td>
<td style="width: 545.933px;">
<p><strong>Description</strong></p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>avg_mu_plot</p>
</td>
<td style="width: 124px;">
<p>Intensity surface</p>
</td>
<td style="width: 545.933px;">
<p>Produces a map of the average intensity for a given subset size</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>compute_intensity</p>
</td>
<td style="width: 124px;">
<p>Intensity surface</p>
</td>
<td style="width: 545.933px;">
<p>Computes the raw and rescaled intensities for a matrix of fitted model coefficients</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>Corr_plot</p>
</td>
<td style="width: 124px;">
<p>Intensity surface</p>
</td>
<td style="width: 545.933px;">
<p>Produces a trace plot of correlation coefficients between the intensity surfaces of the subset models compared to the model fitted with all available points</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>IMSE_plot</p>
</td>
<td style="width: 124px;">
<p>Intensity surface</p>
</td>
<td style="width: 545.933px;">
<p>Produces a trace plot of the integrated mean square error of the intensity surfaces of the subset models compared to the model fitted with all available points</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>makeraster</p>
</td>
<td style="width: 124px;">
<p>Intensity surface</p>
</td>
<td style="width: 545.933px;">
<p>Creates a raster object of a mapped measure from one of the other functions, with an option to export as a .tif file</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>quantilematch</p>
</td>
<td style="width: 124px;">
<p>Intensity surface</p>
</td>
<td style="width: 545.933px;">
<p>Produces a map of misalignment proportions of quantile-categorised intensity surfaces between the subset models of a given subset size and the model fitted with all available points</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>sd_plot</p>
</td>
<td style="width: 124px;">
<p>Intensity surface</p>
</td>
<td style="width: 545.933px;">
<p>Produces a map of standard deviations of the intensity surface for a given subset size</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>coef_plot</p>
</td>
<td style="width: 124px;">
<p>Fitted coefficients</p>
</td>
<td style="width: 545.933px;">
<p>Produces a trace plot of coefficient estimates across models of various subset sizes</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>coef_se_plot</p>
</td>
<td style="width: 124px;">
<p>Fitted coefficients</p>
</td>
<td style="width: 545.933px;">
<p>Produces a trace plot of the standard deviation of the coefficient estimates across models of various subset sizes</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>signcoefs</p>
</td>
<td style="width: 124px;">
<p>Fitted coefficients</p>
</td>
<td style="width: 545.933px;">
<p>Computes the number of positive, zero, and negative coefficient estimates for each covariate across all subset sizes</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>signplot</p>
</td>
<td style="width: 124px;">
<p>Fitted coefficients</p>
</td>
<td style="width: 545.933px;">
<p>Produces a barplot of the estimated coefficient signs for a given covariate across all subset sizes</p>
</td>
</tr>
<tr>
<td style="width: 137px;">
<p>ZeroEnvEffect</p>
</td>
<td style="width: 124px;">
<p>Fitted coefficients</p>
</td>
<td style="width: 545.933px;">
<p>Computes the number of fitted models of each subset size where all coefficients are shrunk to 0</p>
</td>
</tr>
</tbody>
</table>
<p>&nbsp;</p>
