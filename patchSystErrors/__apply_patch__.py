#! python
#vesrion 1.1-3

import os
import sys
import shutil



if sys.version[0] == '2':
    shutil.copy("python2/diff_match_patch.py", "diff_match_patch.py")
if sys.version[0] == '3':
    shutil.copy("python3/diff_match_patch.py", "diff_match_patch.py")

if sys.version[0] == '2':
    subversion = sys.version[2] 
    cmd = '#!/usr/bin/python2.' + str(subversion)+'\n'
    with open("diff_match_patch.py", "r") as infile:
        dmp=infile.readlines()
    dmp[0] = cmd
    outfile = open('diff_match_patch.py','w')
    outfile.writelines(dmp) 
    outfile.close()
   
from diff_match_patch import diff_match_patch, patch_obj

print '''
### based on Diff, Match and Patch Library
###          http://code.google.com/p/google-diff-match-patch/
###          by Neil Fraser
###          Copyright 2006 Google Inc.
'''
print "\n--------------------------"

print "This script will patch your current version of the GSAS-II package \n"
begin = raw_input("Begin [y/n]?") 

if begin=='y':
    fnames = ['GSASIIstrMain.py', 'GSASIIstrMath.py', 'GSASIIobj.py', \
              'GSASIImath.py', 'GSASIIgrid.py', 'GSASII.py', 'config_example.py']   
    fnames0 = ['../'+name for name in fnames]
    fnames1 = ['originalOld/'+name for name in fnames]
    fnames2 = ['modifiedOld/'+name for name in fnames]
    fnames3 = ['originalNew/'+name for name in fnames]
    fnames4 = ['modifiedNew/'+name for name in fnames]
              
    nFiles = len(fnames)

    print "Copying GSASII files to folder 'originalNew':"
    for i in range(nFiles):
        if not os.path.isfile(fnames3[i]):
            print '    ...', fnames[i]
            shutil.copy(fnames0[i], fnames3[i])
        else: 
            print "File", fnames3[i], 'already exists! Please make sure it is the original GSASII file you want to patch. \n'
          

    print "\nGenerating and applying patch:"
    succeed = [False]*nFiles
    for i in range(nFiles):
        diff_obj = diff_match_patch()
        with open (fnames1[i], "r") as infile:
            text1=infile.read()
        with open (fnames2[i], "r") as infile:
            text2=infile.read()     
        with open (fnames3[i], "r") as infile:
            text3=infile.read()    
        pathces = diff_obj.patch_make(text1, text2) 
        text4 = diff_obj.patch_apply(pathces, text3)
        succeed[i] = text4[1]
        text4 = text4[0]
        print "succeeded with",  fnames[i], ":", succeed[i] 
        with open(fnames4[i], "w") as outfile:
            outfile.write(text4)
    
    allSucceeded = True
    for i in range(nFiles):
        if not all(succeed[i]):
            allSucceeded = False
    
    if (allSucceeded):    
        print "\nCopying generated files to GSAS folder:"
        for i in range(nFiles):
            print '    ...', fnames[i]
            shutil.copy(fnames4[i], fnames0[i])
        print "\nPatch was applied successfully!"    
    else:
        print "\nSomething went wrong, sorry"
            
            
end = raw_input("Press RETURN to exit...")