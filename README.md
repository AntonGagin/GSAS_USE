# GSAS_Bayes
Extension to the Rietveld package *[GSAS-II](https://subversion.xor.aps.anl.gov/trac/pyGSAS)* (tested with *GSAS-II* version 0.2.0, revision 1772.) See Toby, B. H. & Von Dreele, R. B. (2013). *J. Appl. Cryst*. **46**, 544-549.

To apply this patch, place **patchSystErrors** folder in *GSAS-II* local folder and run **\_\_apply_patch\_\_**.**py**. If everything is OK, you will see 

```
### based on Diff, Match and Patch Library
###          http://code.google.com/p/google-diff-match-patch/
###          by Neil Fraser
###          Copyright 2006 Google Inc.


--------------------------
This script will patch your current version of GSASII package

Begin [y/n]?
```

Type "y" and follow the instructions. After patch is  completed, start *GSAS-II* normally. In "Controls" menu indicate correction parameters. If several histograms are refined simultaneously, indicate correction parameters, divided by commas, in corresponding order. Set 'E_mu', 'E_beta' or 's' to zero, if you do not want to apply corresponding correction. Do a refinement.
