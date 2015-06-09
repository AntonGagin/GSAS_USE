---
title: "GSAS_USE"
author: "Anton Gagin and Igor Levin"
output: html_document
---

This is an extension to the  *[GSAS-II](https://subversion.xor.aps.anl.gov/trac/pyGSAS)* Rietveld package *GSAS_USE* (Bayesian Statistics Approach to Accounting for <b>U</b>nknown <b>S</b>ystematic <b>E</b>rrors), written and maintained by Anton Gagin (<anton.gagin@nist.gov>.)

*GSAS_USE* addresses the effects of systematic errors in Rietveld refinements. Relevant errors are categorized into multiplicative, additive, and peak-shape types. Corrections for these errors are incorporated into structural refinements using a Bayesian statistics approach, with the corrections themselves treated as nuisance parameters and marginalized out of the analysis. Structural parameters refined using the proposed method represent probability-weighted averages over all possible error corrections. See [Gagin, A. & Levin, I. (2015). *J. Appl. Cryst*. **xx**, xxx-xxx](http://journals.iucr.org/j/) for details.

Current version was tested with *GSAS-II* version 0.2.0, revision 1772.  
For the description of the *GSAS-II* package, see [Toby, B. H. & Von Dreele, R. B. (2013). *J. Appl. Cryst*. **46**, 544-549](http://onlinelibrary.wiley.com/doi/10.1107/S0021889813003531/abstract), or visit their [website](https://subversion.xor.aps.anl.gov/trac/pyGSAS).

***
**Table of Contents**  
**1.** [**Installation**](#install)  
**2.** [**Usage**](#use)  
**3.** [**Description**](#describe)  
**4.** [**Example**](#example)  
**5.** [**What patch does**](#capable)

***

##<a name="install"></a>Installation
To apply this patch, place ***patchSystErrors*** folder in your *GSAS-II* local folder and run **\_\_apply_patch\_\_.py**, or print
```r
execfile(__apply_patch__.py)
```
in python command line interpreter. If everything is OK, you will see 

```r
### based on Diff, Match and Patch Library
###          http://code.google.com/p/google-diff-match-patch/
###          by Neil Fraser
###          Copyright 2006 Google Inc.


--------------------------
This script will patch your current version of the GSAS-II package

Begin [y/n]?
```

Type ```y``` and follow the instructions. 

Folder ***originalOld*** contains some of the original *GSAS-II* source files under revision 1772. Folder ***modifiedOld*** contains our modification of these files. The script copies source files from your current revision of *GSAS-II* into ***originalNew*** folder, so before you apply the patch please make sure that *GSAS-II* local folder contains original *GSAS-II*-files, not our modified versions! **\_\_apply_patch\_\_.py** calculates patch from the difference between files in ***originalOld*** and ***modifiedOld*** folders, applies this patch to the files in ***originalNew*** folder, and writes results to ***modifiedNew*** folder (also to the *GSAS-II* local folder.)

##<a name="use"></a> Usage
After patch is  completed, start *GSAS-II* normally. In **Controls** menu indicate correction parameters. If several histograms are refined simultaneously, indicate correction parameters, divided by commas, in corresponding order. To use similar value for all histograms, put a single number. Set $E\_mu$, $E\_beta$ or $s$ to zero, if you do not want to apply corresponding correction (multiplicative, additive, or peak-shape.) 

If you select *Estimate optimal k\_mu?*, the *Prior factor k\_mu* field will be set to ```optimal```.  The same is true for the *Estimate optimal k\_beta?* and *Prior factor k\_beta* fields.  Deselecting *Estimate optimal k?* will restore the previous value in *Prior factor k*.

If you click on *Correlation length l\_delta* field, the  *estimate it as FWHM /* field will be set to ```none```, and vice versa. The same is also true for the *Stdev sigma\_delta* and  *estimate it as l\_delta/* fields.

To start a Bayesian-corrected refinement, select **Calculate/Refine** in the *GSAS-II* data tree window.

##<a name="describe"></a> Description

* Multiplicative correction $\mu(x)$ is approximated by a set of $E\_mu$ cubic spline functions $\phi_j^{(\mu)}(x)$
$$
\mu(x) = \sum_{j=1}^{E_{\mu}} \left( 1+c_j^{(\mu)}\right) \phi_j^{(\mu)}(x),
$$
where $c_j^{(\mu)}$ are the spline coefficients. Spline knot positions are selected equidistantly.  
Scaling parameter $k\_mu$ indicates how strong is the restriction of closeness of multiplicative correction to unity. It can be estimated from the residual of a standard fit, if *Estimate optimal k\_mu?* is selected.

* Additive correction is approximated by a set of $E\_beta$ cubic spline functions $\phi_j^{(\beta)}(x)$
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
Scaling parameters $sigma\_delta$ and $l\_delta$ indicate standard deviation for the correction and correlation length scale for the point coordinates, respectively. $l\_delta$ can be estimated from the characteristic FWHM values (which depend on x) as $FWHM /p1$, where $p1$ is some number. For a multi-phase problem, if estimated from the FWHM, $l\_delta$ is calculated as a number-of-peaks weighted average over all phases.  
$sigma\_delta$ can be estimated from the $l\_delta$ value(s) as $l\_delta/p2$, where $p2$ is some number.  
To reduce the computation complexity of the problem (you may encounter out-of-memory error for the extremely large histograms) and fasten the calculations, the fitted x-range is divided into $s$ independent segments.
* The iterative procedure works as follows:
	* a standard fit is provided
	* a Bayesian-corrected fit is provided
	* optimal corrections are calculated and applied to the experimental data
	* a Bayesian-corrected fit is provided once again  
As the second Bayesian-corrected fit is provided using the same correction parameters for the already-corrected data, it can easily cause overfitting. Use it only in the case of large systematic errors.

## <a name="example"></a>Example
* [Download](https://subversion.xray.aps.anl.gov/pyGSAS/trunk/help/gsasII.html#Tutorials) exercise files for the 'Combined X-ray/CW-neutron refinement of PbSO4' from *GSAS-II* tutorial. Study the [tutorial](https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/CWCombined/Combined%20refinement.htm) and do a refinement, as described there.
* Deselect all the parameters for refinement except 3 lattice parameters, 11 sets of atomic coordinates, and 5 isotropic atomic displacement parameters. It is especially important that you deselect **Background** and **Histogram scale factor**!
* For this example we want to correct errors of all three types. Set *Number of knots E\_mu* to
```r
15, 20
```
(more splines for the XRD because its residual is worse). These numbers can be further uncreased up to 
```r
30, 45
```
but this will take longer time to calculate. Set *Priot factor k\_mu* to 
```r
1, 1
```
* Set *Number of knots E\_beta* to
```r
15, 20
```
and select *Estimate optimal k\_beta?*
* Set *Number of blocks s* to
```r
8, 8
```
To estimate correlation lengths $l\_delta$ and standard deviations $sigma\_delta$, type
```r
1.5
```
in the *estimate it as FWHM /* and 
```r
2.0
```
in the *estimate it as l\_delta /* fields, respectively.

* Select **Calculate/Refine** in the *GSAS-II* data tree window. The program will provide a standard least-squares fit, followed by the Bayesian-corrected fit. The results will be written, as usual, in the **projectName.lst** file. The details of the Bayesian fit will be stored in the **projectName_cor_iHist.txt** files, where **iHist** is a number indicating histogram index.



## <a name="capable"></a>What patch does

* works with Continuous-Wave and Time-of-Flight data  
* works with single-phase and multi-phase problems  
* calculates optimal values for the prior factors $k\_mu$ and $k\_beta$  
* calculates $FWHM(x)$ dependence to estimate $l\_delta$  
* provides the Bayesian-corrected fits  
* for these fits, calculates most plausible (optimal) corrections  
* stores and draws these corrections   
* draws corresponding residuals and histograms  