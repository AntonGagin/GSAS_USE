import shutil

import os
import sys
import shutil

fnames = ['GSASIIstrMain.py', 'GSASIIstrMath.py', 'GSASIIobj.py', \
          'GSASIImath.py', 'GSASIIgrid.py', 'GSASII.py', 'GSASIIplot.py', \
          'config_example.py']   
fnames0 = ['../'+name for name in fnames]
fnames3 = ['originalNew/'+name for name in fnames]


print "This script will restore original version of the GSAS-II package \n"
begin = raw_input("Begin [y/n]?") 

if begin=='y':
    print "Restoring GSASII files:"
    nFiles = len(fnames)
    for i in range(nFiles):  
        if os.path.isfile(fnames3[i]):
            print '    ...', fnames[i]
            shutil.copy(fnames3[i], fnames0[i])
        else: 
            print "There is no file", fnames[i], ' in \'originalNew\' folder! \n'                         

end = raw_input("\nPress RETURN to exit")