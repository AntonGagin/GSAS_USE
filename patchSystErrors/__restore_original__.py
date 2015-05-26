import shutil




fnames = ['GSASIIstrMain.py', 'GSASIIstrMath.py', 'GSASIIobj.py', \
          'GSASIImath.py', 'GSASIIgrid.py', 'config_example.py']   
fnames0 = ['../'+name for name in fnames]
fnames3 = ['originalNew/'+name for name in fnames]

print "Restoring GSASII files"
nFiles = len(fnames)
for i in range(nFiles):
    print '    ...', fnames[i]
    shutil.copy(fnames3[i], fnames0[i])
    

end = raw_input("Done! Press return to exit")