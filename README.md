# GSAS_Bayes
Extension to the Rietveld package *GSAS-II* (tested with *GSAS-II* version 0.2.0, revision 1772.) See Toby, B. H. & Von Dreele, R. B. (2013). *J. Appl. Cryst*. **46**, 544-549.

Place *__apply_patch__* folder in GSASII local folder and run *__apply_patch__ .py*.

Start GSASII normally. In "Controls" menu indicate correction parameters. If several histograms are refined simultaneously, indicate correction parameters, divided by commas, in corresponding order. Set 'E_mu', 'E_beta' or 's' to zero, if you do not want to apply corresponding correction. Do a refinement.
