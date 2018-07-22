# -*- coding: utf-8 -*-
#config.py - Variables used to set optional configuration options
########### SVN repository information ###################
# $Date: 2018-05-15 20:24:54 +0300 (Tue, 15 May 2018) $
# $Author: vondreele $
# $Revision: 3388 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/config_example.py $
# $Id: config_example.py 3388 2018-05-15 17:24:54Z vondreele $
########### SVN repository information ###################
'''
*config_example.py: Configuration options*
-------------------------------------------

This file contains optional configuration options for GSAS-II. The variables
in this file can be copied to file config.py, which is imported if present.
Access these variables using :func:`GSASIIpath.GetConfigValue`, which returns
None if the variable is not set. Note that a config.py file need not
be present, but if in use it will typically be found with the GSAS-II source
directory (GSASIIpath.Path2GSAS2) or a directory for local GSAS-II
modifications (~/.G2local/ or /Documents and Settings/<User>/.G2local/).  

When defining new config variables for GSAS-II, define them here with a
default value: use None or a string for strings, or use integers or real
values. Include a doc string after each variable is defined to explain
what it does. Use names ending in _location or _directory for items
that will contain directory names.

For example::

    test_int = 0
    test_float = 0.0
    test_string = None (or)
    test_string = 'value'
'''
# </ Anton Gagin
xyFWHM = []  
'''Temporary
Estimated FWHM for Histograms
'''
# Anton Gagin />

debug = False
'''Set to True to turn on debugging mode.This enables use of IPython on
exceptions and on calls to :func:`GSASIIpath.IPyBreak`. Calls to
:func:`GSASIIpath.pdbBreak` will invoke pdb at that location.

If debug is False, calls to :func:`GSASIIpath.IPyBreak` and
:func:`GSASIIpath.pdbBreak` are ignored.
'''

Clip_on = True
''' if True then line plots willl be clipped at plot border; 
if False line plots extend nto white space around plot frme
'''

Transpose = False
'Set to True to cause images to be Transposed when read (for code development)'

Enable_logging = False
'Set to True to enable use of command logging (under development.)'

logging_debug = False
'Set to True to enable debug for logging (under development.)'

Help_mode = "browser"
'''Set to "internal" to use a Python-based web viewer to display
help documentation and tutorials. If set to the default ("browser")
the default web browser is used.
'''

Tutorial_location = None
'''Change this to place tutorials by in a different spot. If None, this defaults to
<user>/My Documents/G2tutorials (on windows) or <user>/G2tutorials. If you want to
use a different location, this can be set here. To install into the location where
GSAS-II is installed, use this::

    Tutorial_location = GSASIIpath.path2GSAS2

As another example, to use ~/.G2tutorials do this::

    Tutorial_location = '~/.G2tutorials'

Note that os.path.expanduser is run on Tutorial_location before it is used.
Also note that GSASIIpath is imported inside config.py; other imports should be
avoided.
'''

Save_paths=False
'''When set to True, the last-used path for saving of .gpx and for
importing of input files is saved in the configuration file.
Note that since this causes the config.py file to be updated whenever files are
saved/imported, any temporary config settings can be saved to disk at that
point.
'''

Starting_directory=None
'''Specifies a default location for starting GSAS-II and where .gpx files
should be read from. Will be updated if Save_paths is True.
Note that os.path.expanduser is run on this before it is used, so the user's
home directory can be specified with a '~'.
'''

Import_directory=None
'''Specifies a default location for importing (reading) input files. Will be
updated if Save_paths is True.
Note that os.path.expanduser is run on this before it is used, so the user's
home directory can be specified with a '~'.
'''

wxInspector = False
'''If set to True, the wxInspector widget is displayed when
GSAS-II is started.
'''

Spot_mask_diameter = 1.0
'''Specifies the default diameter for creation of spot masks. Default is 1.0 mm
'''

Ring_mask_thickness = 0.1
'''Specifies the default thickness for creation of ring and arc masks.
Default is 0.1 degrees 2-theta.
'''

Arc_mask_azimuth = 10.0
'''Specifies the default azimuthal range for creation of arc masks.
Default is 10.0 degrees 2-theta.
'''

Autoint_PollTime = 30.
'''Specifies the frequency, in seconds that AutoInt checks for new files.
Default is 30 seconds
'''

Autoscale_ParmNames = ['userComment2',r'extraInputs\1\extraInputs','Ion_Chamber_I0',]
DefaultAutoScale = "userComment2"
'''Gives the possible selection of incident monitor names as found in an image metadata file.
DefaultAutoScale must be one of the AutoScale_ParmNames
Used in AutoIntegration
'''
Main_Size = '(700,450)'
'''Main window size (width, height) - initially uses wx.DefaultSize but will updated
 and saved as the user changes the window
'''
Main_Pos = '(100,100)'
'''Main window location - will be updated & saved when user moves
it. If position is outside screen then it will be repositioned to default
'''
Plot_Size = '(700,600)'
'''Plot window size (width, height) - initially uses wx.DefaultSize but will updated
 and saved as the user changes the window
'''
Plot_Pos = '(200,200)'
'''Plot window location - will be updated & saved when user moves it
these widows. If position is outside screen then it will be repositioned to default
'''

Tick_length = 8.0
'''Specifies the length of phase tick marks in pixels. Default is 8.'''

Tick_width = 1.0
'''Specifies the width of phase tick marks in pixels.
Fractional values do seem to produce an effect. Default is 1.'''

Contour_color = 'Paired'
''' Specifies the color map to be used for contour plots (images, pole figures, etc.)
will be applied for new images and if Saved for a new start of GSAS-II
'''

Movie_fps = 10
''' Specifies movie frames-per-second; larger number will make smoother modulation movies but larger files.
'''

Movie_time = 5
''' Specifices time in sec for one modulation loop; larger number will give more frames for same fps'
'''

Multiprocessing_cores = 0
''' Specifies the number of cores to use when performing multicore computing. A number less
than zero causes the recommended number of cores [using multiprocessing.cpu_count()/2]
to be used. Setting this number to 0 or 1 avoids use of the multiprocessing module: all
computations are performed in-line. 
'''

Show_timing = False
'''If True, shows various timing results.'''

Column_Metadata_directory = None
'''When specified and when images are read, GSAS-II will read metadata from a 1-ID
style .par and a .EXT_lbls (EXT = image extension) or .lbls file. See :func:`GSASIIIO.readColMetadata` for
information on how this is done.
'''

Instprm_default = False
'''when True, GSAS-II instprm file are shown as default; when False, old GSAS stype prm, etc files are default
'''

Plot_Colors = 'k r g b m c'
'''The colors for line plots: use one of 'k'-black, 'r'-red, 'b'-blue, 'g'-green, 'm'-magenta, 'c'-cyan for the
line colors in order of obs., calc., back., diff., color5 & color6 separated by spaces; 6 items required.
'''