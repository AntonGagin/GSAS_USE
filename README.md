# GSAS_Bayes
Extension to the Rietveld package *[GSAS-II](https://subversion.xor.aps.anl.gov/trac/pyGSAS)* (tested with *GSAS-II* version 0.2.0, revision 1772.) See Toby, B. H. & Von Dreele, R. B. (2013). *J. Appl. Cryst*. **46**, 544-549.

## Installation
To apply this patch, place **patchSystErrors** folder in *GSAS-II* local folder and run **\_\_apply_patch\_\_.py**. If everything is OK, you will see 

```
### based on Diff, Match and Patch Library
###          http://code.google.com/p/google-diff-match-patch/
###          by Neil Fraser
###          Copyright 2006 Google Inc.


--------------------------
This script will patch your current version of the GSAS-II package

Begin [y/n]?
```

Type "y" and follow the instructions. 

Folder **originalOld** contains some of the original *GSAS-II* source files under revision 1772. Folder **modifiedOld** contains our modification of these files. The script copies source files from your current revision of *GSAS-II* into **originalNew** folder, so before you apply the patch please make sure that *GSAS-II* local folder contains original *GSAS-II*-files, not our modified versions! **\_\_apply_patch\_\_.py** calculates patch from the difference between files in **originalOld** and **modifiedOld** folders, applies this patch to the files in **originalNew** folder, and writes results to **modifiedNew** folder (also to the *GSAS-II* local folder.)

## Usage
After patch is  completed, start *GSAS-II* normally. In "Controls" menu indicate correction parameters. If several histograms are refined simultaneously, indicate correction parameters, divided by commas, in corresponding order. Set $E\_mu$, $E\_beta$ or $s$ to zero, if you do not want to apply corresponding correction. Do a refinement.

## Description
GSAS_Bayes addresses the effects of systematic errors in Rietveld refinements. Relevant errors are categorized into multiplicative, additive, and peak-shape types. Corrections for these errors are incorporated into structural refinements using a Bayesian statistics approach, with the corrections themselves treated as nuisance parameters and marginalized out of the analysis. Structural parameters refined using the proposed method represent probability-weighted averages over all possible error corrections. 

* Multiplicative correction is approximated by a set of $E\_mu$ cubic spline functions $\phi_j^{(\mu)}(x)$: 
$$
\mu(x) = \sum_{j=1}^{E_{\mu}} \left( 1+c_j^{(\mu)}\right) \phi_j^{(\mu)}(x).
$$
Scaling parameter $k\_mu$ indicates how strong is the restriction of closeness of multiplicative correction to unity.

* Additive correction is approximated by a set of $E\_beta$ cubic spline functions $\phi_j^{(\beta)}(x)$: 
$$
\beta(x) = \sum_{j=1}^{E_{\beta}} c_j^{(\beta)} \phi_j^{(\beta)}(x).
$$
Scaling parameter $k\_beta$ indicates how strong is the smoothness restriction on the additive correction.

* A diffraction profile is corrected by moving points of a diffraction curve in the horizontal direction. A probability of each 'move' $\delta x$ is calculated as 
$$
p(\delta x) \propto \exp \left(  -\frac{1}{2}\delta x^T \Sigma_{\delta}^{-1} \delta x \right),
$$
where the covariance matrix $\Sigma_{\delta}^{-1}$ is defined as
$$
\Sigma_{ij}^{(\delta)} = \sigma_{\delta}^2 \exp \left(  -\frac{1}{2} \left( \frac{x_i-x_j}{l_{\delta}} \right)^2 \right).
$$
Scaling parameters $sigma\_delta$ and $l\_delta$ indicate standard deviation for the correction and correlation length scale for the point coordinates, respectively. These can be estimated from the characteristic FWHM values (which depend on x).
To reduce the computation complexity of the problem (you may encounter out-of-memory error for the extremely large histograms) and fasten the calculations, the fitted x-range is divided into $s$ independent segments.

## Example
* [Download](https://subversion.xray.aps.anl.gov/pyGSAS/trunk/help/gsasII.html#Tutorials) exercise files for the 'Combined X-ray/CW-neutron refinement of PbSO4' from *GSAS-II* tutorial. Study the [tutorial](https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/CWCombined/Combined%20refinement.htm) and do a refinement, as described there.
* Deselect all the parameters for refinement except 3 lattice parameters, 11 sets of atomic coordinates, and 5 isotropic atomic displacement parameters. It is especially important that you deselect Background and Histogram scale factor!
* For this example we want to correct errors of all three types. Set 'Number of knots E\_mu' to
```r
15, 20
```
(more splines for the XRD because its residual is worse). These numbers can be further uncreased up to 
```r
30, 45
```
but this will take longer time to calculate. Set 'Priot factor k\_mu' to 
```r
1, 1
```
* Set 'Number of knots E\_beta' to
```r
15, 20
```
and select 'Estimate optimal k\_beta?'
* Set 'Number of blocks s' to
```r
8, 8
```
To estimate correlation lengths $l\_delta$ and standard deviations $sigma\_delta$, type
```r
1.5
```
in the 'estimate it as FWHM /' and 
```r
2.0
```
in the 'estimate it as FWHM /' fields, respectively.
* Select **Calculate/Refine** in the GSAS-II data tree window. The program will provide a standard least-squares fit, followed by the Bayesian-corrected fit. The results will be written, as usual, in the *project\_name.lst* file. The details of the Bayesian fit will be stored in the *project\_name_cor_iHist.txt* files, where iHist indicates histogram index.



