import shutil

fnames = ['GSASIIstrMain.py', 'GSASIIstrMath.py', 'GSASIIobj.py', \
          'GSASIImath.py', 'GSASIIgrid.py', 'config_example.py']   
fnames0 = ['../'+name for name in fnames]
fnames3 = ['originalNew/'+name for name in fnames]

print "Restoring GSASII files:"
nFiles = len(fnames)
for i in range(nFiles):  
    if os.path.isfile(fnames3[i]):
        print '    ...', fnames[i]
        shutil.copy(fnames3[i], fnames0[i])
    else: 
        print "There is no file", fnames[i], ' in \'originalNew\' folder! \n'
                       
end = raw_input("\nDone! Press RETURN to exit")