---
title: "GSAS_USE"
author: "Anton Gagin and Igor Levin"
output: html_document
---

**Please open README.pdf to see the formulas**

This is an extension to the  *[GSAS-II](https://subversion.xor.aps.anl.gov/trac/pyGSAS)* Rietveld package *GSAS_USE* (Bayesian Statistics Approach to Accounting for <b>U</b>nknown <b>S</b>ystematic <b>E</b>rrors), written and maintained by Anton Gagin (av.gagin@gmail.com, igor.levin@nist.gov) 

*GSAS_USE* addresses the effects of systematic errors in Rietveld refinements. The errors are categorized into multiplicative, additive, and peak-shape types. Corrections for these errors are incorporated into using a Bayesian statistics approach, with the corrections themselves treated as nuisance parameters and marginalized out of the analysis. Structural parameters refined using the proposed method represent probability-weighted averages over all possible error corrections. See [Gagin, A. & Levin, I. (2015). *Accounting for Unknown Systematic Errors in Rietveld Refinements: A Bayesian Statistics Approach.* *J. Appl. Cryst*. **48**, 1201-1211](http://journals.iucr.org/j/issues/2015/04/00/po5042/stdsup.html) for details.

The current version has been tested with *GSAS-II* version 0.2.0, revision 1970.  
For details of the *GSAS-II* package, refer to [Toby, B. H. & Von Dreele, R. B. (2013). *J. Appl. Cryst*. **46**, 544-549](http://onlinelibrary.wiley.com/doi/10.1107/S0021889813003531/abstract), or visit their [website](https://subversion.xor.aps.anl.gov/trac/pyGSAS).

***
**Table of Contents**  
**1.** [**Installation**](#installation)  
**2.** [**Usage**](#usage)  
**3.** [**Description**](#description)  
**4.** [**Example**](#example)  
**5.** [**Bugs**](#bugs)  

***

##<a name="installation"></a>Installation

To apply this patch, place the ***patchSystErrors*** folder in your *GSAS-II* local folder and run **\_\_apply_patch\_\_.py**, or print
```r
execfile('__apply_patch__.py')
```
in a python command line interpreter. If everything works correctly, the following message will be displayed 

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

Folder ***originalOld*** contains some of the original *GSAS-II* source files under revision 1970. Folder ***modifiedOld*** contains our modification of these files. The script copies the source files from your current revision of *GSAS-II* into the ***originalNew*** folder. Before applying the patch please ensure that the local folder with *GSAS-II* contains the original *GSAS-II*-files and not the modified versions! **\_\_apply\_patch\_\_.py** calculates the patch from the difference between the files in the ***originalOld*** and ***modifiedOld*** folders, applies this patch to the files in the ***originalNew*** folder, and writes the results to the ***modifiedNew*** folder (as well as to the *GSAS-II* local folder.) 

To restore the original *GSAS-II*-files, run **\_\_restore\_original\_\_.py**.

To update patch, run **\_\_update\_patch\_\_.py**.


##<a name="usage"></a> Usage
After the patch has been  applied, start *GSAS-II* normally. In **Controls** menu specify the correction parameters. If several histograms are refined simultaneously, list these parameters, separated by commas, in the order corresponding to the order of the histograms (it may not correspond to their order on the data tree). If you wish to the same value of the parameter for all histograms, enter a single number. Set $E\_mu$, $E\_beta$ or $s$ to zero, if you do not want to apply a particular correction (multiplicative, additive, or peak-shape.) 

If you select *Estimate optimal k\_mu?*, the *Prior factor k\_mu* field will be set to ```optimal```.  The same is true for the *Estimate optimal k\_beta?* and *Prior factor k\_beta* fields.  Deselecting *Estimate optimal k?* will restore the previous value in *Prior factor k*.

If you click on *Correlation length l\_delta* field, the  *estimate it as FWHM /* field will be set to ```none```, and vice versa. The same is true for the fields *Stdev sigma\_delta* and  *estimate it as l\_delta/*.

To start a Bayesian-corrected refinement, select **Calculate/Refine** in the *GSAS-II* data tree window. To see refinement results, select **Data/Open .lst file** or **Data/Compare standard and Bayesian fits**.

##<a name="description"></a> Description

* The multiplicative correction $\mu(x)$ is approximated by a set of $E_{\mu}$ cubic spline functions $\phi_j^{(\mu)}(x)$
$$
\mu(x) = \sum_{j=1}^{E_{\mu}} \left( 1+c_j^{(\mu)}\right) \phi_j^{(\mu)}(x),
$$
where $c_j^{(\mu)}$ are the spline coefficients. Spline-knot positions are selected equidistantly.  
The scaling parameter $k_{\mu}$ reflects the strength of the restriction on closeness of the multiplicative correction to unity. It can be estimated by the program from the residual of a standard fit (no corrections), if *Estimate optimal k\_mu?* is selected.

* The additive correction is approximated using a set of $E_{\beta}$ cubic spline functions $\phi_j^{(\beta)}(x)$
$$
\beta(x) = \sum_{j=1}^{E_{\beta}} c_j^{(\beta)} \phi_j^{(\beta)}(x).
$$
The scaling parameter $k_{\beta}$ reflects the strength of the smoothness restriction on the additive correction.

* A diffraction profile is corrected by varying x-coordinates of the individual points of a diffraction curve. A probability of each 'move' $\delta x$ is calculated as 
$$
p(\delta x) \propto \exp \left(  -\frac{1}{2}\delta x^T \Sigma_{\delta}^{-1} \delta x \right),
$$
where the covariance matrix $\Sigma_{\delta}^{-1}$ is defined as
$$
\Sigma_{ij}^{(\delta)} = \sigma_{\delta}^2 \exp \left(  -\frac{1}{2} \left( \frac{x_i-x_j}{l_{\delta}} \right)^2 \right).
$$
The scaling parameters $\sigma_{\delta}$ and $l_{\delta}$ describe a standard deviation for the correction and correlation length for the point coordinates, respectively. $l_{\delta}$ can be estimated from characteristic FWHM values for diffraction peaks (which depend on x) as $FWHM /p1$, where $p1$ can be any real number. For a multi-phase refinement, if estimated from FWHM, $l_{\delta}$ is calculated as a an average weighted by a number of peaks for all the phases. Fig. 1 provides a hint on how to select $p1$ for $l_{\delta}$.
 
<div style="width:450px; height=450px">
![Figure 1](https://cloud.githubusercontent.com/assets/8290742/9686784/97321530-52f3-11e5-9a7b-adf22a7b24f8.png)
</div>

$\sigma_{\delta}$ can be estimated from the $l_{\delta}$ value(s) as $l_{\delta}/p2$, where $p2$ can be any real number. Normally, $p2 \approx 1.5-2$
To reduce computational complexity (e.g. one may get an out-of-memory error for extremely large histograms) and speed the calculations up, the fitted x-range is divided into $s$ independent segments. 

* The iterative procedure works as follows:
  
    * a standard fit is performed
    * a Bayesian-corrected fit is performed
    * the optimal corrections are calculated and applied to the experimental data
    * a Bayesian-corrected fit is repeated
  
The second Bayesian-corrected fit is prone to overfitting because it uses the same correction parameters as those that have been already applied to the data. Therefore, we advise to limit the use of the iterative option to cases of large systematic errors.

* If your select *run sampler for MCMC?* the patch will do the following:
	* perform a standard fit
	* call the [*emcee*](http://dan.iel.fm/emcee/current/) library and run the Goodman & Weare's Affine Invariant MCMC sampler
	* perform a Bayesian-corrected fit to obtain the final estimates
  
Results of the MCMC sampler will be saved in a text file and as a plot in a project folder. Prior to using this feature make sure that [*emcee*](http://dan.iel.fm/emcee/current/user/install/) and [*triangle_plot*](https://github.com/dfm/triangle.py) libraries are installed.
  
## <a name="example"></a>Example
* [Download](https://subversion.xray.aps.anl.gov/pyGSAS/trunk/help/gsasII.html#Tutorials) the example files for a 'Combined X-ray/CW-neutron refinement of PbSO4' from the *GSAS-II* tutorial. Perform the refinements as described in the [tutorial](https://subversion.xray.aps.anl.gov/pyGSAS/Tutorials/CWCombined/Combined%20refinement.htm).
* Deselect all the refinable parameters except for the structural variables which include 3 lattice parameters, 11 sets of atomic coordinates, and 5 isotropic atomic displacement parameters. MAKE SURE to deselect **Background** and **Histogram scale factor**!
* For this example we want to correct all three types of errors. Set the *Number of knots E\_mu* to
```r
15, 20
```
(more splines are selected for the XRD data because it exhibits the worse residual). These numbers of knots can be increased up to 
```r
30, 45
```
but this will take longer to calculate. Set *Prior factor k\_mu* to 
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

* Select **Calculate/Refine** in the *GSAS-II* data tree window. The program will perform a standard least-squares fit followed by a Bayesian-corrected fit. The results will be saved in the **projectName.lst** file. The details of the Bayesian fit will be stored in the **projectName_cor_iHist.txt** files, where **iHist** is the histogram number.

Select **Data/Open .lst file** to see the *GSAS-II* .lst project file. The residuals are summarized in the table entitled as

```r
********************************************************
*
* == SUMMARIZING REFINEMENT RESULTS: ==
```

Calculated as a sum of squared residuals for the Bayesian approach are expected to be larger than those obtained using standard LS technique. Calculated with optimal corrections residuals are expected to be smaller. 

Select **Data/Compare standard and Bayesian fits** to see fit results. The notation for the parameters is the following:
```r
i::Name:j
```
Here $i$ and $j$ indicate histogram and atom number, respectively, and $Name$ indicates parameter name. Note, that *GSAS-II* fits the changes in atomic coordinates rather than their absolute values. These changes are calculated with respect to the starting values. Absolute values for the atomic coordinates are given in the .lst project file.

## <a name="bugs"></a>Bugs
To report a bug or ask a question, send an e-mail to both of us (<av.gagin@gmail.com> and <igor.levin@nist.gov>). For a bug report, please include the error message and traceback from the console window [text beginning with "Traceback (most recent call..."].

Please cite Gagin, A. & Levin, I. (2015). *Accounting for Unknown Systematic Errors in Rietveld Refinements: A Bayesian Statistics Approach.* *J. Appl. Cryst*. **48**, 1201-1211 in publications that use this method.

