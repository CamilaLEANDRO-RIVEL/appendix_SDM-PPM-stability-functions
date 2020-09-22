# appendix_SDM-PPM-stability-functions
This repository contains the R code to (i) run a species distribution model in the point process modelling framework and to (ii) explore model stability (or the capacity to predict correctly all presence data) within the PPM perspective fitted with a lasso penalty and observer bias corrections. Different R functions which were written to establish intensity and coefficient measures. These functions are the stability assessment toolbox, referred to as “diagnostic tools”. Code, simulated data and a detailed tutorial illustrating use of this code are provided. The model is described in a paper submitted to Ecological Modelling.

Hereafter: A short description of the supplied R functions to explore model stability, which are contained in the DiagnosticFunctions.R file supplied in the supplementary material. Full details of these functions and a demonstration of their usage appears in the RMarkdown tutorial in the supplementary material. 
 
 <table>
<tbody>
<tr>
<td width="137">
<p>Function</p>
</td>
<td width="124">
<p>Characteristic</p>
</td>
<td width="305">
<p>Description</p>
</td>
</tr>
<tr>
<td width="137">
<p>avg_mu_plot</p>
</td>
<td width="124">
<p>Intensity surface</p>
</td>
<td width="305">
<p>Produces a map of the average intensity for a given subset size</p>
</td>
</tr>
<tr>
<td width="137">
<p>compute_intensity</p>
</td>
<td width="124">
<p>Intensity surface</p>
</td>
<td width="305">
<p>Computes the raw and rescaled intensities for a matrix of fitted model coefficients</p>
</td>
</tr>
<tr>
<td width="137">
<p>Corr_plot</p>
</td>
<td width="124">
<p>Intensity surface</p>
</td>
<td width="305">
<p>Produces a trace plot of correlation coefficients between the intensity surfaces of the subset models compared to the model fitted with all available points</p>
</td>
</tr>
<tr>
<td width="137">
<p>IMSE_plot</p>
</td>
<td width="124">
<p>Intensity surface</p>
</td>
<td width="305">
<p>Produces a trace plot of the integrated mean square error of the intensity surfaces of the subset models compared to the model fitted with all available points</p>
</td>
</tr>
<tr>
<td width="137">
<p>makeraster</p>
</td>
<td width="124">
<p>Intensity surface</p>
</td>
<td width="305">
<p>Creates a raster object of a mapped measure from one of the other functions, with an option to export as a .tif file</p>
</td>
</tr>
<tr>
<td width="137">
<p>quantilematch</p>
</td>
<td width="124">
<p>Intensity surface</p>
</td>
<td width="305">
<p>Produces a map of misalignment proportions of quantile-categorised intensity surfaces between the subset models of a given subset size and the model fitted with all available points</p>
</td>
</tr>
<tr>
<td width="137">
<p>sd_plot</p>
</td>
<td width="124">
<p>Intensity surface</p>
</td>
<td width="305">
<p>Produces a map of standard deviations of the intensity surface for a given subset size</p>
</td>
</tr>
<tr>
<td width="137">
<p>coef_plot</p>
</td>
<td width="124">
<p>Fitted coefficients</p>
</td>
<td width="305">
<p>Produces a trace plot of coefficient estimates across models of various subset sizes</p>
</td>
</tr>
<tr>
<td width="137">
<p>coef_se_plot</p>
</td>
<td width="124">
<p>Fitted coefficients</p>
</td>
<td width="305">
<p>Produces a trace plot of the standard deviation of the coefficient estimates across models of various subset sizes</p>
</td>
</tr>
<tr>
<td width="137">
<p>signcoefs</p>
</td>
<td width="124">
<p>Fitted coefficients</p>
</td>
<td width="305">
<p>Computes the number of positive, zero, and negative coefficient estimates for each covariate across all subset sizes</p>
</td>
</tr>
<tr>
<td width="137">
<p>signplot</p>
</td>
<td width="124">
<p>Fitted coefficients</p>
</td>
<td width="305">
<p>Produces a barplot of the estimated coefficient signs for a given covariate across all subset sizes</p>
</td>
</tr>
<tr>
<td width="137">
<p>ZeroEnvEffect</p>
</td>
<td width="124">
<p>Fitted coefficients</p>
</td>
<td width="305">
<p>Computes the number of fitted models of each subset size where all coefficients are shrunk to 0</p>
</td>
</tr>
</tbody>
</table>
<p>&nbsp;</p>
