import shutil

import os
import sys
import shutil
import urllib
import zipfile

f = open('__apply_patch__.py', 'r')
version = f.readline()
version = f.readline()

url = 'https://github.com/AntonGagin/GSAS_USE/archive/master.zip'
filehandle, _ = urllib.urlretrieve(url)
zf = zipfile.ZipFile(filehandle, 'r')
data = zf.read("GSAS_USE-master/patchSystErrors/__apply_patch__.py")

i1 = data.find('\n', 0) + 1
i2 = data.find('\n', i1) + 1
if(data[i1:i2]==version):
    print 'Your patch is already up to date! \n'
    end = raw_input("\nPress RETURN to exit")
    sys.exit()
    
    
print "You are going to do the following steps: \n 1. Restore original GSAS-II file \n 2. Download patch updates \n 3. Apply patch updates \n"
begin = raw_input("Begin [y/n]?") 

if begin=='y': 
    print 'Step 1: \n'
    execfile("__restore_original__.py")
    
    print 'Step 2: \n'
    print 'This script will update your GSAS_USE patch \n'
    begin2 = raw_input("Begin [y/n]?") 
    if begin2=='y':
        url = 'https://github.com/AntonGagin/GSAS_USE/archive/master.zip'
        filehandle, _ = urllib.urlretrieve(url)
        zf = zipfile.ZipFile(filehandle, 'r')

        for filename in zf.namelist():
            if ("GSAS_USE-master/patchSystErrors/__apply_patch__.py" in filename):      
                data = zf.read(filename)
                name = filename.replace('GSAS_USE-master/patchSystErrors/','')
                print '   updating ', name, '...'
                f = open(name, 'w')
                f.write(data)
                f.close()      
            if ("GSAS_USE-master/patchSystErrors/__restore_original__.py" in filename):      
                data = zf.read(filename)
                name = filename.replace('GSAS_USE-master/patchSystErrors/','')
                print '   updating ', name, '...'
                f = open(name, 'w')
                f.write(data)
                f.close()          
            if ("GSAS_USE-master/patchSystErrors/originalOld/" in filename):
                if(".py" in filename):
                    data = zf.read(filename)
                    name = filename.replace('GSAS_USE-master/patchSystErrors/originalOld/','')
                    name = os.path.join('originalOld/', name)
                    print '   updating ', name, '...'
                    f = open(name, 'w')
                    f.write(data)
                    f.close()
            if ("GSAS_USE-master/patchSystErrors/modifiedOld/" in filename):
                if(".py" in filename):
                    data = zf.read(filename)
                    name = filename.replace('GSAS_USE-master/patchSystErrors/modifiedOld/','')
                    name = os.path.join('modifiedOld/', name)
                    print '   updating ', name, '...'
                    f = open(name, 'w')
                    f.write(data)
                    f.close()

    print ' \n'                
        
    print 'Step 3: \n'
    execfile("__apply_patch__.py")

   
end = raw_input("\nPress RETURN to exit")
