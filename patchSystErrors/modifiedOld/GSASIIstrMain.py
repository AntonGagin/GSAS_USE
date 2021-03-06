# -*- coding: utf-8 -*-
'''
*GSASIIstrMain: main structure routine*
---------------------------------------

'''
########### SVN repository information ###################
# $Date: 2018-07-13 22:44:01 +0300 (Fri, 13 Jul 2018) $
# $Author: toby $
# $Revision: 3471 $
# $URL: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASIIstrMain.py $
# $Id: GSASIIstrMain.py 3471 2018-07-13 19:44:01Z toby $
########### SVN repository information ###################
from __future__ import division, print_function
import platform
import sys
import os.path as ospath
import time
import math
import copy
if '2' in platform.python_version_tuple()[0]:
    import cPickle
else:
    import pickle as cPickle
import numpy as np
import numpy.linalg as nl
import scipy.optimize as so
import GSASIIpath
GSASIIpath.SetBinaryPath()
GSASIIpath.SetVersionNumber("$Revision: 3471 $")
import GSASIIlattice as G2lat
import GSASIIspc as G2spc
import GSASIImapvars as G2mv
import GSASIImath as G2mth
import GSASIIstrIO as G2stIO
import GSASIIstrMath as G2stMth
import GSASIIobj as G2obj

sind = lambda x: np.sin(x*np.pi/180.)
cosd = lambda x: np.cos(x*np.pi/180.)
tand = lambda x: np.tan(x*np.pi/180.)
asind = lambda x: 180.*np.arcsin(x)/np.pi
acosd = lambda x: 180.*np.arccos(x)/np.pi
atan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi

ateln2 = 8.0*math.log(2.0)
DEBUG = True

# </ Anton Gagin
#import scipy as sp
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import config_example                
from scipy import interpolate
from scipy import linalg
from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy.ma as ma


def RefineCore(Controls,Histograms,Phases,restraintDict,rigidbodyDict,parmDict,varyList,
    calcControls,pawleyLookup,ifSeq,printFile,dlg,GPXfile):
    '''Core optimization routines, shared between SeqRefine and Refine

    :returns: 5-tuple of ifOk (bool), Rvals (dict), result, covMatrix, sig
    '''
#    print 'current',varyList
#    for item in parmDict: print item,parmDict[item] ######### show dict just before refinement
    G2mv.Map2Dict(parmDict,varyList)
    ifPrint = True
    if ifSeq:
        ifPrint = False
    Rvals = {}
    histoList = Histograms.keys()
    histoList.sort()
    
# we'll need it later, but this cycle has to be done before the first refinement
    nHist = 0
    hIdList = []
    for histogram in histoList:
        Histogram = Histograms[histogram]
        hId = Histogram['hId']
        hIdList.append(hId) 
        nHist = nHist + 1
    config_example.xyFWHM = [[0]*(max(hIdList)+1), [0]*(max(hIdList)+1)]
    
# Anton Gagin />	
    while True:
        begin = time.time()
        values =  np.array(G2stMth.Dict2Values(parmDict, varyList))
        if np.any(np.isnan(values)):
            raise G2obj.G2Exception('ERROR - nan found in LS parameters - use Calculate/View LS parms to locate')
        # test code to compute GOF and save for external repeat
        #args = ([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg)
        #print '*** before fit chi**2',np.sum(G2stMth.errRefine(values,*args)**2)
        #fl = open('beforeFit.cpickle','wb')
        #cPickle.dump(values,fl,1)
        #cPickle.dump(args[:-1],fl,1)
        #fl.close()
        Ftol = Controls['min dM/M']
        Xtol = Controls['SVDtol']
        Factor = Controls['shift factor']
        if 'Jacobian' in Controls['deriv type']:
            result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = int(result[2]['nfev']/2)
        elif 'analytic Hessian' in Controls['deriv type']:
            Lamda = Controls.get('Marquardt',-3)
            maxCyc = Controls['max cyc']
            result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,xtol=Xtol,maxcyc=maxCyc,Print=ifPrint,lamda=Lamda,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = result[2]['num cyc']+1
            Rvals['lamMax'] = result[2]['lamMax']
            Controls['Marquardt'] = -3  #reset to default
        elif 'Hessian SVD' in Controls['deriv type']:
            maxCyc = Controls['max cyc']
            result = G2mth.HessianSVD(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,xtol=Xtol,maxcyc=maxCyc,Print=ifPrint,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = result[2]['num cyc']+1
        else:           #'numeric'
            result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
            ncyc = 1
            if len(varyList):
                ncyc = int(result[2]['nfev']/len(varyList))
#        table = dict(zip(varyList,zip(values,result[0],(result[0]-values))))
#        for item in table: print item,table[item]               #useful debug - are things shifting?
        runtime = time.time()-begin
        Rvals['SVD0'] = result[2].get('SVD0',0)
        Rvals['converged'] = result[2].get('Converged')
        Rvals['DelChi2'] = result[2].get('DelChi2',-1.)
        Rvals['chisq'] = np.sum(result[2]['fvec']**2)
        G2stMth.Values2Dict(parmDict, varyList, result[0])
        G2mv.Dict2Map(parmDict,varyList)
        Rvals['Nobs'] = Histograms['Nobs']
        Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histograms['sumwYo'])*100.      #to %
        Rvals['GOF'] = np.sqrt(Rvals['chisq']/(Histograms['Nobs']-len(varyList)))
        printFile.write(' Number of function calls: %d No. of observations: %d No. of parameters: %d User rejected: %d Sp. gp. extinct: %d\n'%  \
            (result[2]['nfev'],Histograms['Nobs'],len(varyList),Histograms['Nrej'],Histograms['Next']))
        printFile.write(' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles\n'%(runtime,runtime/ncyc,ncyc))
        printFile.write(' wR = %7.2f%%, chi**2 = %12.6g, GOF = %6.2f\n'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF']))
        sig = len(varyList)*[None,]
        if 'None' in str(type(result[1])) and ifSeq:    #this bails out of a sequential refinement on singular matrix
            IfOK = False
            covMatrix = []
            print ('**** Refinement failed - singular matrix ****')
            if 'Hessian' in Controls['deriv type']:
                num = len(varyList)-1
                for i,val in enumerate(np.flipud(result[2]['psing'])):
                    if val:
                        print ('Bad parameter: '+varyList[num-i])
            else:
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        print ('Bad parameter: '+varyList[ipvt-1])
            break
        IfOK = True
        try:
            covMatrix = result[1]*Rvals['GOF']**2
            sig = np.sqrt(np.diag(covMatrix))
            if np.any(np.isnan(sig)) or not sig.shape:
                print ('*** Least squares aborted - some invalid esds possible ***')
#            table = dict(zip(varyList,zip(values,result[0],(result[0]-values)/sig)))
#            for item in table: print item,table[item]               #useful debug - are things shifting?
            break                   #refinement succeeded - finish up!
        except TypeError:          #result[1] is None on singular matrix or LinAlgError
            IfOK = False
            if not len(varyList):
                covMatrix = []
                break
            print ('**** Refinement failed - singular matrix ****')
            if 'Hessian' in Controls['deriv type']:
                num = len(varyList)-1
                for i,val in enumerate(np.flipud(result[2]['psing'])):
                    if val:
                        print ('Removing parameter: '+varyList[num-i])
                        del(varyList[num-i])
            else:
                Ipvt = result[2]['ipvt']
                for i,ipvt in enumerate(Ipvt):
                    if not np.sum(result[2]['fjac'],axis=1)[i]:
                        print ('Removing parameter: '+varyList[ipvt-1])
                        del(varyList[ipvt-1])
                        break
						
# </ Anton Gagin
    val_GSAS = result[0]
    std_GSAS = sig       
    GOF_GSAS = Rvals['GOF']
    wR_GSAS = Rvals['Rwp']
# =1. INITIALIZATION=
# =1.1 Controls=        
# E_mu: number of splines for multiplicative factor
# optK_mu: find optimal k_mu
# k_mu: prior parameter for multiplicative factor
    E_mu = Controls['corrParam E_mu'].split(',')
    try:
        E_mu = [int(p) for p in E_mu]*nHist
        if 'EstimateKMu' in Controls:
            optK_mu = Controls['EstimateKMu']
        else: 
            optK_mu = False
        k_mu = [0]*nHist   
        if(not optK_mu):
            k_mu = Controls['corrParam k_mu'].split(',')
            k_mu = [float(p) for p in k_mu]*nHist
    except (ValueError, RuntimeError, TypeError, NameError):
        E_mu = [0]*nHist  	
# E_beta: number of splines for additive factor
# optK_beta: find optimal k_beta
# k_beta: prior parameter for additive factor
    E_beta = Controls['corrParam E_beta'].split(',')
    try:
        E_beta = [int(p) for p in E_beta]*nHist
        if 'EstimateKBeta' in Controls:
            optK_beta = Controls['EstimateKBeta']
        else: 
            optK_beta = False
        k_beta = [0]*nHist   
        if(not optK_beta):
            k_beta = Controls['corrParam k_beta'].split(',')
            k_beta = [float(p) for p in k_beta]*nHist
    except (ValueError, RuntimeError, TypeError, NameError):
        E_beta = [0]*nHist  		
# sigma_delta: standard deviation for peak-shape correction
# l_delta: correlation length for peak-shape correction
    l_delta = Controls['corrParam l_delta'].split(',')
    try:
        xL_delta = []
        if l_delta[0] == 'none':
            divN = Controls['corrParam FWHMDivN'].split(',')   # l_delta = fwhm/divN
            divN = [float(p) for p in divN]*nHist
            divN = np.array(divN)[range(nHist)]  
            xL_delta = [config_example.xyFWHM[0][p] for p in hIdList]
            l_delta = [config_example.xyFWHM[1][p] for p in hIdList]        
            l_delta = [np.array(l_delta[i])/divN[i] for i in range(nHist)]     
        else:        
            l_delta = [np.array(float(p)) for p in l_delta]*nHist 
    except (ValueError, RuntimeError, TypeError, NameError):
        l_delta = [0]*nHist  	
    
    sigma_delta = Controls['corrParam sigma_delta'].split(',')
    try:	
        if sigma_delta[0] == 'none':
            divN = Controls['corrParam l_deltaDivN'].split(',')
            divN = [float(p) for p in divN]*nHist
            divN = np.array(divN)[range(nHist)] 
            sigma_delta = [np.array(l_delta[i])/divN[i] for i in range(nHist)]
        else:        
            sigma_delta = [np.array(float(p)) for p in sigma_delta]*nHist
    except (ValueError, RuntimeError, TypeError, NameError):
        sigma_delta = [0]*nHist  	
# nBlocks: number of blocks s to use for each Histogram 
    nBlocks = Controls['corrParam num blocks s'].split(',')
    try:
        nBlocks = [int(p) for p in nBlocks]*nHist
    except (ValueError, RuntimeError, TypeError, NameError):
        nBlocks = [0]*nHist  	
# Markov chain Monte Carlo
    doMCMC = Controls['doMCMC']
    
    if(np.any(E_beta) or np.any(E_mu) or np.any(nBlocks)):
        print ("Bayesian-corrected refinement")
# =1.2 Functions=
#
        def SE(xi, xj, sd, l):
            res = sd**2*np.exp(-0.5*(xi-xj)**2/l**2)
            return(res)     
#            
        def basisMatrix(x, knots_x, der=0):
            E = len(knots_x)
            N = len(x)
            bM = np.zeros(shape=(N,E))
            for i in range(1,E+1,1):
                bM[:, i-1] = basisSpline(x, knots_x, i, der)
            return(bM)
#            
        def basisSpline(x, knots_x, knots_i, der=0):
            knots_y = 0*np.array(knots_x)
            knots_y[knots_i-1] = 1.0
            tck = interpolate.splrep(knots_x, knots_y)
            y = interpolate.splev(x, tck, der)
            return(y)
#            
        def getSpline(x, knots_x, knots_y, der=0):
            bM = basisMatrix(x, knots_x)
            spl = np.dot(bM,knots_y)
            return(spl)
#            
        def overlap(x, knots_x, knots_i, knots_j, der=0):
            res = np.multiply(basisSpline(x, knots_x, knots_i, der), basisSpline(x, knots_x, knots_j, der))
            return res   
#            
        def blockMult(K, nBlocks, sBlock, b):     
            b = np.array(b)
            x0 = 0
            x1 = 0
            if (K == []):
                res = b
            else: 
                if (b.ndim==1):        
                    res=[]
                    for iBlock in range(nBlocks):
                        x1 += sBlock[iBlock]
                        res = np.concatenate((res, np.dot(K[iBlock], b[x0:x1])))       # K*b             
                        x0 += sBlock[iBlock]
                else:        
                    m = b.shape[1]        
                    for iBlock in range(nBlocks):
                        x1 += sBlock[iBlock]
                        if (iBlock==0):
                            res = np.dot(K[iBlock], b[x0:x1, range(m)]) 
                        else:
                            res = np.concatenate((res, np.dot(K[iBlock], b[x0:x1, range(m)])))               
                        x0 += sBlock[iBlock] 
            return(res)                       
#           
        def extrap1d(interpolator):
            xs = interpolator.x
            ys = interpolator.y
            def pointwise(x):
                if x < xs[0]:
                    return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
                elif x > xs[-1]:
                    return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
                else:
                    return interpolator(x)
            def ufunclike(xs):
                return np.array(map(pointwise, np.array(xs)))
            return ufunclike  
# =1.3 Additional objects=  
#   
# sBlock: size of each block for each histogram
# histLen: xF-xB
# yp: y derivative
# mDE_mu: overlap matrix of second derivatives for multiplicative factor
# mDE_beta: overlap matrix for additive factor
# mPhiC, mPhiB: Phi-matrices for multiplicative and additive factors, respectively
# mInvMC: inverse of mMC
# mMC: M-matrix for multiplicative correction     
# mInvMB: inverse of mMB
# mMB: M-matrix for additive correction
# mK: resulting (main) covariance matrix
# mK_temp: I - mK

        sBlock =   [0]*nHist
        histLen =  [0]*nHist
        yp =       [0]*nHist
        mDE_mu =     [0]*nHist
        mDE_beta =     [0]*nHist
        mPhiC =    [0]*nHist
        mPhiB =    [0]*nHist
        mInvMC =   [np.zeros(shape=(0, 0)) for i in range(nHist)]
        mMC =      [np.zeros(shape=(0, 0)) for i in range(nHist)]
        mInvMB =   [np.zeros(shape=(0, 0)) for i in range(nHist)]
        mMB =      [np.zeros(shape=(0, 0)) for i in range(nHist)]
        mK =      [[] for i in range(nHist)]   # mK = [ [], [], ... ]
        mK_temp = [0]*nHist
        hId = -1
# =2. PREPARING CORRECTIONS=
        for histogram in histoList:
            hId = hId + 1
            Histogram = Histograms[histogram]
            hIdG = Histogram['hId']
            hfx = ':%d:'%(hIdG)
            wtFactor = calcControls[hfx+'wtFactor']
            Limits = calcControls[hfx+'Limits']
            x,y,w,yc,yb,yd = Histogram['Data']
            yc *= 0.0                           #zero full calcd profiles
            yb *= 0.0
            yd *= 0.0
            xB = np.searchsorted(x,Limits[0])
            xF = np.searchsorted(x,Limits[1])+1
            yc[xB:xF],yb[xB:xF] = G2stMth.getPowderProfile(parmDict,x[xB:xF],
                varyList,Histogram,Phases,calcControls,pawleyLookup)
            yc[xB:xF] += yb[xB:xF]
            yd[xB:xF] = y[xB:xF]-yc[xB:xF]
            weights = list(wtFactor*w)
            weights = np.sqrt(weights[xB:xF])
            histLen[hId] = xF-xB

            mPhiC[hId] = np.zeros(shape=(histLen[hId], E_mu[hId]))
            mPhiB[hId] = np.zeros(shape=(histLen[hId], E_beta[hId]))
                    
    #####
    # SIGMA^-1 = mK - mK*Phi*mInvMC*Phi.T*mK  // mult + shape error together
    # SIGMA^-1 = I - Phi*mInvMC*Phi.T         // multiplicative error only
    # SIGMA^-1 = mK                           // peak shape error only

#####        
# =2.1 Marginalization over multiplicative error=
            print ("\n=== Histogram #", hId+1, "===")
            if (E_mu[hId]):  
                x_min = x[xB]
                x_max = x[xF-1]
                knots_x = x_min + (x_max-x_min)/(E_mu[hId]-1)*np.array(range(E_mu[hId]))            
                mDE_mu[hId] = np.zeros(shape=(E_mu[hId],E_mu[hId]))
                for knots_i in range(E_mu[hId]):
                    for knots_j in range(E_mu[hId]):
                        mDE_mu[hId][knots_i, knots_j] = quad(overlap, x_min, x_max, args=(knots_x, knots_i, knots_j))[0]
                mPhiC[hId] = basisMatrix(x[xB:xF], knots_x)     
    #            gQ = list(yc[xB:xF]-yb[xB:xF])
                gQ = list(yc[xB:xF])
                for i in range(mPhiC[hId].shape[0]):
                    mPhiC[hId][i,:] = mPhiC[hId][i,:]*weights[i]*gQ[i]
                
    # Find optimal k_mu </       
                if(optK_mu):
                    print ("Calculating optimal k_mu:")
                    f = np.multiply(yd[xB:xF], weights)
                    c_opt = 0.2*np.random.rand(E_mu[hId])-0.1  # initial guess
                    c_opt_tmp = c_opt
                    ccd = np.array([0.005]*E_mu[hId])         # to avoid null argument in log
                    delta = np.dot(ccd.T, np.dot(mDE_mu[hId], ccd))
#                    print ("delta", delta)
                    itermax = 20
                    for i in range(itermax):
                        MM = np.dot(mPhiC[hId].T, mPhiC[hId]) + E_mu[hId]*mDE_mu[hId]/(np.dot((c_opt).T, np.dot(mDE_mu[hId], c_opt)) + delta)
                        b = np.dot(mPhiC[hId].T, f)        						
                        c_opt = linalg.solve(ma.getdata(MM), ma.getdata(b))
                        print ("  iteration #", i, "c_opt", c_opt)
                        if ( (np.mean(abs(c_opt_tmp-c_opt)) < 1e-4) and (i > 3)):
                            print ("  converged!")
                            break
                        c_opt_tmp = c_opt
                        
                    cc_opt = np.array(getSpline(x[xB:xF], knots_x, c_opt))
                    dx = x[xB+1]-x[xB]
                    dx = np.diff(x[xB:xF])
                    dx = np.append(dx,dx[-1])
                    bb = np.array(cc_opt*cc_opt)
                    integr = sum(bb*dx)
                    k_mu[hId] = 2.0/integr
                    print ("\n  optimal k_mu=", str(k_mu[hId]))
    # /> Find optimal k_mu            
                # MM not yet inverse -- see below

#####
#  =2.2 Marginalization over multiplicative error=
            if (E_beta[hId]):  
                x_min = x[xB]
                x_max = x[xF-1]
                knots_x = x_min + (x_max-x_min)/(E_beta[hId]-1)*np.array(range(E_beta[hId]))            
                mDE_beta[hId] = np.zeros(shape=(E_beta[hId],E_beta[hId]))
                der=2
                for knots_i in range(1, E_beta[hId]+1):
                    for knots_j in range(1,E_beta[hId]+1):
                        mDE_beta[hId][knots_i-1, knots_j-1] = quad(overlap, x_min, x_max, args=(knots_x, knots_i, knots_j, der))[0]
               
                mPhiB[hId] = basisMatrix(x[xB:xF], knots_x)     
                for i in range(mPhiB[hId].shape[0]):
                    mPhiB[hId][i,:] = mPhiB[hId][i,:]*weights[i]
                
    # Find optimal k_beta        
                if(optK_beta):
                    print ("Calculating optimal k_beta:")
                    f = np.multiply(yd[xB:xF], weights)
                    b_opt = 0.2*np.random.rand(E_beta[hId])-0.1  # initial guess
                    b_opt_temp = b_opt  # initial guess
                    ccd = np.array([0.005]*E_beta[hId])         # to avoid null argument in log
                    delta = np.dot(ccd.T, np.dot(mDE_beta[hId], ccd))
                    
                    itermax = 20
                    for i in range(itermax):
                        MM = np.dot(mPhiB[hId].T, mPhiB[hId]) + E_beta[hId]*mDE_beta[hId]/(np.dot((b_opt).T, np.dot(mDE_beta[hId], b_opt)) + delta)
                        b = np.dot(mPhiB[hId].T, f)
                        b_opt = linalg.solve(ma.getdata(MM), ma.getdata(b))
                        print ("  iteration #", i, "b_opt", b_opt)
                        if ( (np.mean(abs(b_opt_temp-b_opt)) < 1) and (i > 3)):
                            print ("  converged!")
                            break
                        b_opt_temp = b_opt
                        
                        
                    bb_opt = np.array(getSpline(x[xB:xF], knots_x, b_opt, der))
                    dx = np.diff(x[xB:xF])
                    dx = np.append(dx,dx[-1])
                    bb = np.array(bb_opt*bb_opt)
                    integr = sum(bb*dx)
                    k_beta[hId] = 2.0/integr
                    print ("\n  optimal k_beta=", str(k_beta[hId]))
    # / Find optimal k_beta        
                
#####
#  =2.3 Marginalization over peak-shape error=
            if(not nBlocks[hId]):
                sBlock[hId] = [histLen[hId]]   
            else:
                print ("Calculating peak-shape correction:")
                yp[hId] = np.zeros(shape=histLen[hId])  # y prime
                dx = np.diff(x[xB:xF])
                dx = np.append(dx,dx[-1])
                for i in range(1, histLen[hId]):
                    yp[hId][i]=(yc[i+xB]-yc[i-1+xB])/dx[i]
                yp[hId][0]=yp[hId][1]
                yp[hId] = np.multiply(yp[hId], weights)

                mK[hId]=[0]*nBlocks[hId]
                
                sBlock[hId] = np.array([0]*nBlocks[hId])
                sBlock[hId][range(nBlocks[hId])] = histLen[hId]/nBlocks[hId]
                sBlock[hId][nBlocks[hId]-1] += histLen[hId] % nBlocks[hId]
                       
                for iBlock in range(nBlocks[hId]):
                    mK[hId][iBlock]=np.zeros(shape=(sBlock[hId][iBlock],sBlock[hId][iBlock]))
                                 
                print ("  Calculating K...")
                x0 = xB   
                if(xL_delta==[]):
                    for iBlock in range(nBlocks[hId]):
                        print ("    step %s of %s" % (iBlock+1, nBlocks[hId]))
                        ld2 = -0.5/l_delta[hId]**2.              
                        sd2 = sigma_delta[hId]**2.                    
                        for i in range(sBlock[hId][iBlock]):
                            for j in range(sBlock[hId][iBlock]):
                                mK[hId][iBlock][i,j] = sd2*np.exp(ld2*(x[i+x0]-x[j+x0])**2)*yp[hId][i+x0-xB]*yp[hId][j+x0-xB]                           
                        x0 += sBlock[hId][iBlock]                
                else:
                    if (sigma_delta[hId].shape==()):
                        sigma_delta[hId] = [sigma_delta[hId]]*len(xL_delta[hId])    
                    
                    xL_delta[hId], l_delta[hId], sigma_delta[hId] = \
                        (list(t) for t in zip(*sorted(zip(xL_delta[hId], l_delta[hId], sigma_delta[hId]))))

                    xL_delta[hId] = np.array(xL_delta[hId])
                    l_delta[hId] = np.array(l_delta[hId])
                    sigma_delta[hId] = np.array(sigma_delta[hId])            
                    
                    nanInd = np.where(np.isnan(l_delta[hId]))
                    xL_delta[hId] = np.delete(xL_delta[hId], nanInd)
                    l_delta[hId] = np.delete(l_delta[hId], nanInd)
                    sigma_delta[hId] = np.delete(sigma_delta[hId], nanInd)
                    
                    ld_func = interp1d(xL_delta[hId], l_delta[hId])
                    ld_func = extrap1d(ld_func)
                    ld = ld_func(x[xB:xF])  

                    ld = np.sqrt(2.)*np.array(ld)
                    ld[np.isnan(ld)] = np.mean(l_delta[hId])                    
                                        
                    
                    p = np.polyfit(x[xB:xF], ld, 3)
                    p = np.poly1d(p)
                    p[p < 1e-6] = np.min(l_delta[hId])
                    ld = p(x[xB:xF])
                    ld[ld < 1e-6] = np.min(l_delta[hId])            
                    ld2 = [0.5*p**2. for p in ld]
                     
                    sd_func = interp1d(xL_delta[hId], sigma_delta[hId])
                    sd_func = extrap1d(sd_func)
                    sd = sd_func(x[xB:xF])     
                    sd[np.isnan(sd)] = np.mean(sigma_delta[hId])
                    
                    p = np.polyfit(x[xB:xF], sd, 3)
                    p = np.poly1d(p)
                    p[p < 1e-6] = np.min(sigma_delta[hId])
                    sd = p(x[xB:xF])     
                    sd[sd < 1e-6] = np.min(sigma_delta[hId])                  
                    sd2 = [p**2. for p in sd]   
                    
                    for iBlock in range(nBlocks[hId]):
                        print ("    step %s of %s" % (iBlock+1, nBlocks[hId]))
                        for i in range(sBlock[hId][iBlock]):
                            for j in range(sBlock[hId][iBlock]):
                                mK[hId][iBlock][i,j] = 0.5*(sd2[i+x0-xB] + sd2[j+x0-xB]) \
                                                       *np.sqrt(ld[i+x0-xB]*ld[j+x0-xB]/(ld2[i+x0-xB]+ld2[j+x0-xB])) \
                                                       *np.exp( -(x[i+x0]-x[j+x0])**2./(ld2[i+x0-xB]+ld2[j+x0-xB]))  \
                                                       *yp[hId][i+x0-xB]*yp[hId][j+x0-xB]             
                        x0 += sBlock[hId][iBlock]
                           
                mK_temp[hId]=[0]*nBlocks[hId]
                for iBlock in range(nBlocks[hId]):
                    mK_temp[hId][iBlock]=mK[hId][iBlock] + np.identity(sBlock[hId][iBlock])
                # we'll need mK_temp later; mK_temp = K*(K+I)**(-1)
                #                           mK = I - K*(K+I)**(-1)
                
                print ("  Calculating (K+I)**(-1)...")
                for iBlock in range(nBlocks[hId]):         
                    mK_temp[hId][iBlock] = nl.inv(mK_temp[hId][iBlock])
                
                print ("  Calculating K*(K+I)**(-1)...")
                for iBlock in range(nBlocks[hId]):         
                    mK_temp[hId][iBlock] = np.dot(mK[hId][iBlock], mK_temp[hId][iBlock])     
                for iBlock in range(nBlocks[hId]):         
                    mK[hId][iBlock] = np.identity(sBlock[hId][iBlock]) - mK_temp[hId][iBlock]
#####
#  =2.4 Combining corrections=                 
            print ("Calculating M-matrix...")
            if E_mu[hId]:
                # M mult
                print ("  multiplicative")
                mMC[hId] = np.dot(mPhiC[hId].T, mPhiC[hId]) + 2*k_mu[hId]*mDE_mu[hId]    
                # M mult + peak shape
                if nBlocks[hId]:            
                    print ("  multiplicative + peak-shape")
                    mMC[hId] = mMC[hId] - np.dot(mPhiC[hId].T, blockMult(mK_temp[hId], nBlocks[hId], sBlock[hId], mPhiC[hId]))
                mInvMC[hId] = nl.inv(mMC[hId])

            if E_beta[hId]:
                # M add
                print ("  additive")
                mInvMB[hId] = np.dot(mPhiB[hId].T, mPhiB[hId]) + 2*k_beta[hId]*mDE_beta[hId]  
                # M add + peak shape
                if nBlocks[hId]:            
                    print ("  additive + peak-shape")
                    mInvMB[hId] = mInvMB[hId] - np.dot(mPhiB[hId].T, blockMult(mK_temp[hId], nBlocks[hId], sBlock[hId], mPhiB[hId]))
                mMB[hId]=mInvMB[hId]    
                # M add + mult
                if E_mu[hId]:    
                    print ("  multiplicative + additive")
                    matr_1 = np.dot(mPhiC[hId].T, mPhiB[hId])  # PhiC.T*PhiB 
                    matr_1 = np.dot(mInvMC[hId], matr_1)  # MC^-1*PhiC.T*PhiB      
                    matr_1 = np.dot(mPhiC[hId], matr_1)   # PhiC*MC^-1*PhiC.T*PhiB      
                    mInvMB[hId] = mInvMB[hId] - np.dot(mPhiB[hId].T, matr_1)
                # M add + mult + peak-shape
                if nBlocks[hId] & E_mu[hId]:
                    print ("  multiplicative + additive + peak-shape")
                    matr_2 =   blockMult(mK_temp[hId], nBlocks[hId], sBlock[hId], matr_1)  
                    matr_2 = np.dot(mPhiB[hId].T,matr_2)
                    matr_3 =   blockMult(mK_temp[hId], nBlocks[hId], sBlock[hId], mPhiB[hId])              
                    matr_3 = np.dot(mPhiC[hId].T, matr_3) 
                    matr_3 = np.dot(mInvMC[hId], matr_3) 
                    matr_3 = np.dot(mPhiC[hId], matr_3) 
                    matr_3 =   blockMult(mK_temp[hId], nBlocks[hId], sBlock[hId], matr_3)              
                    matr_3 = np.dot(mPhiB[hId].T, matr_3) 
                    mInvMB[hId] = mInvMB[hId] + matr_2 + matr_2.T - matr_3
                mInvMB[hId] = nl.inv(mInvMB[hId])    
                    
          
        #data['marg mu'] = ', '.join(map(str,mu))        
        iHist=0
        mPhiCTotal = mPhiC[0]
        mPhiBTotal = mPhiB[0]
        mInvMCTotal = mInvMC[0]
        mInvMBTotal = mInvMB[0]
#####
#  =2.5 Joint matrices for all histograms=         
        print ("Calculating joint matrices...")
        while (iHist < (nHist-1)):  
            # mPhiC1 =| mPhiC1    0     0   ... |
            #         |  0      mPhiC2  0   ... |
            #         |  0        0  mPhiC3 ... |
            #         |  0       ...            |   
            n1 = mPhiC[iHist+1].shape[0]
            n2 = mPhiCTotal.shape[1]
            P1 = np.concatenate((mPhiCTotal,np.zeros([n1, n2])))
            
            n1 = mPhiCTotal.shape[0]
            n2 = mPhiC[iHist+1].shape[1]            
            P2 = np.concatenate((np.zeros([n1, n2]), mPhiC[iHist+1]))
            mPhiCTotal = np.concatenate((P1.T, P2.T)).T   
            
            # mPhiB1 =| mPhiB1    0       0   ... |
            #         |  0      mPhiB2    0   ... |
            #         |  0        0    mPhiB3 ... |
            #         |  0          ...           |   
            n1 = mPhiB[iHist+1].shape[0]
            n2 = mPhiBTotal.shape[1]
            P1 = np.concatenate((mPhiBTotal,np.zeros([n1, n2])))
            
            n1 = mPhiBTotal.shape[0]
            n2 = mPhiB[iHist+1].shape[1]            
            P2 = np.concatenate((np.zeros([n1, n2]), mPhiB[iHist+1]))
            mPhiBTotal = np.concatenate((P1.T, P2.T)).T           
                 
            # MC1 = | MC1   0  ...|
            #       |  0   MC2 ...|
            #       |     ...     |     
            n1 = mInvMC[iHist+1].shape[0]
            n2 = mInvMCTotal.shape[1]        
            I1 = np.concatenate((mInvMCTotal, np.zeros([n1, n2])))
            
            n1 = mInvMCTotal.shape[0]
            n2 = mInvMC[iHist+1].shape[1]         
            I2 = np.concatenate((np.zeros([n1, n2]), mInvMC[iHist+1]))
            mInvMCTotal = np.concatenate((I1.T, I2.T)).T   

            # MB1 = | MB1   0  ...|
            #       |  0   MB2 ...|
            #       |     ...     |     
            n1 = mInvMB[iHist+1].shape[0]
            n2 = mInvMBTotal.shape[1]        
            I1 = np.concatenate((mInvMBTotal, np.zeros([n1, n2])))
            
            n1 = mInvMBTotal.shape[0]
            n2 = mInvMB[iHist+1].shape[1]         
            I2 = np.concatenate((np.zeros([n1, n2]), mInvMB[iHist+1]))
            mInvMBTotal = np.concatenate((I1.T, I2.T)).T           
            
            iHist = iHist + 1  

#####
#  =3 REFINEMENT!=              
        peakCor={'doCor':False}
        if np.any(nBlocks):
            peakCor = {'doCor':True, 'nBlocks':nBlocks, 'sBlock':sBlock, 'nHist':nHist, 'mK':mK} 
        multCor={'doCor':False}
        if np.any(E_mu):
            multCor={'doCor':True, 'mPhiC':mPhiCTotal, 'mInvMC':mInvMCTotal}
        addCor={'doCor':False}
        if np.any(E_beta):
            addCor={'doCor':True, 'mPhiB':mPhiBTotal, 'mInvMB':mInvMBTotal}      
            
        while True:            
####                
# =3.0 MCMC=
            if doMCMC:
                print ('Running MCMC...')
                nwalkers = int(Controls['nwalkers'])
                nIterMCMC = int(Controls['nIterMCMC'])
                if (nwalkers==0): nwalkers=50
                if (nIterMCMC==0): nIterMCMC=100
                sampler = G2mth.MarginalizedMCMC(G2stMth.errRefine,values, 2*sig,
                                    peakCor=peakCor, multCor=multCor, addCor=addCor,
                                    args=([Histograms,Phases,restraintDict,rigidbodyDict],
                                    parmDict,varyList,calcControls,pawleyLookup,dlg), 
                                    nwalkers=nwalkers, nIter=nIterMCMC, Print=False)
                ndim=len(sig)     
                nToOmit = min(nIterMCMC/2, 50)                
                samples = sampler.chain[:, nToOmit:, :].reshape((-1, ndim))                                          
                bins = min(nIterMCMC*nwalkers/5, 20)
                fnametxt = ospath.splitext(GPXfile)[0]+'-MCMC.txt'
                np.savetxt(fnametxt, sampler.flatchain, header=str(varyList))
                
                try:
                    import triangle
                    fig = triangle.corner(samples, bins=bins, labels=varyList, 
                                          quantiles=[0.16, 0.5, 0.84])
                    fnamepng = ospath.splitext(GPXfile)[0]+'-MCMC.png'
                    fig.savefig(fnamepng)
                except:
                    print ('Cannot allocate memory for the image. Image saved as', fnamepng)
                
                if 'Jacobian' in Controls['deriv type']:            
                    result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                        ftol=Ftol,col_deriv=True,factor=Factor,
                        args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                elif 'Hessian' in Controls['deriv type']:
                    result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,maxcyc=maxCyc,Print=ifPrint,
                        args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
                else:           #'numeric'
                    result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                        args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
    
                
# =3.1 Bayesian fit =
            maxCyc = Controls['max cyc']
            result = G2mth.MarginalizedLSQ(G2stMth.errRefine,values,Jac=G2stMth.dervRefine,  
                peakCor=peakCor, multCor=multCor, addCor=addCor, optCor={},
                ftol=Ftol,maxcyc=maxCyc, Print=True, args=([Histograms,Phases,restraintDict,rigidbodyDict],
                parmDict,varyList,calcControls,pawleyLookup,dlg))      
                
            G2stMth.Values2Dict(parmDict, varyList, result[0])
            G2mv.Dict2Map(parmDict,varyList)
     
            (x, b_opt, bb_opt, c_opt, cc_opt, dx_opt, cb_opt,
                ydiff, ystd, yexp, yuncor, ycor) =  ([0]*nHist for i in range(12))
            f_opt  = []                        
                    
            print ("Calculating optimal corrections")
            hId = -1 
            for histogram in histoList:
                hId = hId + 1
                Histogram = Histograms[histogram]
                hIdG = Histogram['hId']
                hfx = ':%d:'%(hIdG)
                wtFactor = calcControls[hfx+'wtFactor']
                Limits = calcControls[hfx+'Limits']
                x[hId],y,w,yc,yb,yd = Histogram['Data']
                yc *= 0.0                           #zero full calcd profiles
                yb *= 0.0
                yd *= 0.0
                xB = np.searchsorted(x[hId],Limits[0])
                xF = np.searchsorted(x[hId],Limits[1])+1
                x[hId] = x[hId][xB:xF]
                yc[xB:xF],yb[xB:xF] = G2stMth.getPowderProfile(parmDict,x[hId],
                    varyList,Histogram,Phases,calcControls,pawleyLookup)
                yc[xB:xF] += yb[xB:xF]
                yd[xB:xF] = y[xB:xF]-yc[xB:xF]
                weights = list(wtFactor*w)
                weights = np.sqrt(weights[xB:xF])
                f = np.multiply(yd[xB:xF], weights)
# =3.2 Calculate optimal corrections=    
                c_opt[hId] = np.zeros(shape=(0,))
                cc_opt[hId] = [1.]*histLen[hId]  # null correction is 1
                b_opt[hId] = np.zeros(shape=(0,))
                bb_opt[hId] = [0.]*histLen[hId]  # null correction is 0
                cb_opt[hId] = np.zeros(shape=(0,))
                if(E_mu[hId] or E_beta[hId]):
                    cb_opt[hId] = blockMult(mK[hId], nBlocks[hId], sBlock[hId], f)
                    cb_opt[hId] = np.dot(np.concatenate((mPhiC[hId].T, mPhiB[hId].T)), cb_opt[hId])  # len = E_beta+E_mu
                    if(E_mu[hId] and E_beta[hId]):
                        mMCB = np.dot(mPhiC[hId].T, blockMult(mK[hId], nBlocks[hId], sBlock[hId], mPhiB[hId]))
                    else:
                        mMCB = np.zeros(shape=(mMC[hId].shape[0], mMB[hId].shape[1]))
                    mMtot1 = np.concatenate((mMC[hId].T, mMCB.T)).T
                    mMtot2 = np.concatenate((mMCB, mMB[hId].T)).T
                    mMtotInv = nl.inv(np.concatenate((mMtot1, mMtot2)))
                    cb_opt[hId] = np.dot(mMtotInv, cb_opt[hId])
                    
                    x_min = x[hId][0]
                    x_max = x[hId][-1]
                    if(E_mu[hId]):
                        c_opt[hId] = cb_opt[hId][0:E_mu[hId]] + 1.
                        knots_x = x_min + (x_max-x_min)/(E_mu[hId]-1)*np.array(range(E_mu[hId]))
                        cc_opt[hId] = np.array(getSpline(x[hId], knots_x, c_opt[hId]))
                    if(E_beta[hId]):
                        b_opt[hId] = cb_opt[hId][E_mu[hId]:(E_beta[hId]+E_mu[hId])]
                        knots_x = x_min + (x_max-x_min)/(E_beta[hId]-1)*np.array(range(E_beta[hId]))
                        bb_opt[hId] = np.array(getSpline(x[hId], knots_x, b_opt[hId]))                    
                      
                dx_opt[hId] = [0]*histLen[hId]
                if(nBlocks[hId]):
                    y2 = np.dot(mPhiC[hId], c_opt[hId]-1.) + np.dot(mPhiB[hId], b_opt[hId])
                    y2 = f - y2
                    dx_opt[hId] = blockMult(mK_temp[hId], nBlocks[hId], sBlock[hId], y2)
                    dx_opt[hId] = np.multiply(1.0/yp[hId], dx_opt[hId])   
         
# =3.3 Apply optimal corrections=    
                yuncor[hId] = np.zeros(shape=(xF-xB))
                yuncor[hId][0:(xF-xB)] = yc[xB:xF]            
                x_cor = x[hId]-dx_opt[hId]
                yc_func = interp1d(x_cor, yc[xB:xF])
                yc_func = extrap1d(yc_func)
                yc[xB:xF] = yc_func(x[hId])
                
   #            yc[xB:xF] = np.multiply(cc_opt[hId], (yc[xB:xF] - yb[xB:xF])) + yb[xB:xF]
                yc[xB:xF] = np.multiply(cc_opt[hId], yc[xB:xF])
                yc[xB:xF] = yc[xB:xF] + bb_opt[hId]
                yd[xB:xF] = y[xB:xF]-yc[xB:xF]
                
                ydiff[hId] = yd[xB:xF]
                ystd[hId] = 1./weights
                ycor[hId] = yc[xB:xF]
                yexp[hId] = y[xB:xF]
                f_opt = np.concatenate((f_opt, -np.multiply(yd[xB:xF], weights)))
                        
            GOF_uncor = np.sqrt(np.sum(result[2]['fvec']**2)/(Histograms['Nobs']-len(varyList)))  
            wR_uncor = np.sqrt(np.sum(result[2]['fvec']**2)/Histograms['sumwYo'])*100.
            optCor = {'dx_opt':dx_opt, 'cc_opt':cc_opt, 'bb_opt':bb_opt}
            
#####
# =4 ITERATION=
            doIter = Controls['doIter']
            if doIter:
                print ('Second Iteration:')
                result = G2mth.MarginalizedLSQ(G2stMth.errRefine_opt,values,Jac=G2stMth.dervRefine,  
                    peakCor=peakCor, multCor=multCor, addCor=addCor, optCor=optCor,
                    ftol=Ftol, maxcyc=maxCyc, Print=True, args=([Histograms,Phases,restraintDict,rigidbodyDict],
                    parmDict,varyList,calcControls,pawleyLookup,dlg))         
                G2stMth.Values2Dict(parmDict, varyList, result[0])
                G2mv.Dict2Map(parmDict,varyList)
                

                
####                
# =5 SAVING RESULTS=
            result[2]['fvec'] = f_opt
            result[2]['chisq1'] = np.sum(f_opt**2)    
            ncyc = result[2]['num cyc']+1
            Rvals['lamMax'] = result[2]['lamMax'] 
            runtime = time.time()-begin
            Rvals['converged'] = result[2].get('Converged')
            Rvals['DelChi2'] = result[2].get('DelChi2',-1.)
            Rvals['chisq'] = np.sum(result[2]['fvec']**2)
            Rvals['Nobs'] = Histograms['Nobs']
            Rvals['Rwp'] = np.sqrt(Rvals['chisq']/Histograms['sumwYo'])*100.      #to %
            Rvals['GOF'] = np.sqrt(Rvals['chisq']/(Histograms['Nobs']-len(varyList)))
            GOF_opt = Rvals['GOF']
            wR_opt = Rvals['Rwp']

            print ('GOF standard approach: ', GOF_GSAS)
            print ('GOF before corrections: ', GOF_uncor)
            print ('GOF after after corrections: ', GOF_opt)
            
            printFile.write(' Applying corrections: \n')
            line = ' Number of function calls:' + str(result[2]['nfev']) + ' Number of observations: ' + str(Histograms['Nobs']) + ' Number of parameters: ' + str(len(varyList)) + '\n'
            printFile.write(line)
            printFile.write(' Refinement time = %8.3fs, %8.3fs/cycle, for %d cycles\n'%(runtime,runtime/ncyc,ncyc))
            printFile.write(' wR = %7.2f%%, chi**2 = %12.6g, GOF = %6.2f\n'%(Rvals['Rwp'],Rvals['chisq'],Rvals['GOF']))
            
            printFile.write( '\n'  )
            printFile.write(90*'*')
            printFile.write( '*\n')
            printFile.write( '* == SUMMARIZING REFINEMENT RESULTS: ==\n')
            printFile.write( '*\n')
            printFile.write( '* standard least-squares approach:        wR = %7.2f%%, reduced chi**2 = %6.2f\n'%(wR_GSAS,GOF_GSAS**2))
            printFile.write( '* Bayesian approach:                      wR = %7.2f%%, reduced chi**2 = %6.2f\n'%(wR_uncor,GOF_uncor**2))
            printFile.write( '* Bayesian approach, optimal corrections: wR = %7.2f%%, reduced chi**2 = %6.2f\n'%(wR_opt,GOF_opt**2))
            printFile.write( '*\n')
            printFile.write(90*'*')
            printFile.write( '\n'  )
            
            
            IfOK = True
            try:
                covMatrix = result[1]*GOF_GSAS**2
                sig = np.sqrt(np.diag(covMatrix))
                if np.any(np.isnan(sig)):
                    print ('*** Least squares aborted - some invalid esds possible ***')
                break                   #refinement succeeded - finish up!
            except TypeError,FloatingPointError:          #result[1] is None on singular matrix
                IfOK = False
                if not len(varyList):
                    covMatrix = []
                    sig = []
                    break
                print ('**** Refinement failed - singular matrix ****')
                if 'Hessian' in Controls['deriv type']:
                    num = len(varyList)-1
                    for i,val in enumerate(np.flipud(result[2]['psing'])):
                        if val:
                            print ('Removing parameter: ',varyList[num-i])
                            del(varyList[num-i])                    
                else:
                    Ipvt = result[2]['ipvt']
                    for i,ipvt in enumerate(Ipvt):
                        if not np.sum(result[2]['fjac'],axis=1)[i]:
                            print ('Removing parameter: ',varyList[ipvt-1])
                            del(varyList[ipvt-1])
                            break
      
        histNames=G2stIO.GetHistogramNames(GPXfile,['PWDR',])    
        from multiprocessing import Process
        p = Process(target=plotCorrections, args=(nHist, histNames, E_mu, E_beta, nBlocks, x, cc_opt, bb_opt, dx_opt, ydiff, ystd, yexp, ycor, yuncor))
        p.start()        
              
      
        val_Bayes = result[0]
        std_Bayes = sig
        rownames = np.array(varyList)
        printFileCompare = open(ospath.splitext(GPXfile)[0]+'-compare.lst','w')
        header = '''    Name   val(GSAS) sig(GSAS) val(Bayes) sig(Bayes)'''       
        dat = np.array([val_GSAS, std_GSAS, val_Bayes, std_Bayes  ]).T
        printFileCompare.write(header)
        for i in range(len(rownames)):
            printFileCompare.write('%10s % 1.2e % 1.2e % 1.2e % 1.2e' % (rownames[i], dat[i,0], dat[i,1], dat[i,2], dat[i,3]))       
        printFileCompare.close()


        for i in range(nHist):
            header = histNames[i]+'''\n        x        mult         add      dtheta        resid       stdev        yexp        ycor      yuncor'''
            printCorFile = ospath.splitext(GPXfile)[0]+'_cor_'+str(i)+'.txt'
            dat = np.array([x[i], cc_opt[i], bb_opt[i], dx_opt[i], ydiff[i], ystd[i], yexp[i], ycor[i], yuncor[i]])        
            np.savetxt(fname=printCorFile, X=dat.T, header=header, fmt='% 1.5e')
    
    elif(doMCMC):
        print ('Running MCMC...')
        nwalkers = int(Controls['nwalkers'])
        nIterMCMC = int(Controls['nIterMCMC'])
        if (nwalkers==0): nwalkers=50
        if (nIterMCMC==0): nIterMCMC=100
        sampler = G2mth.StandardMCMC(G2stMth.errRefine,values, 2*sig,
                            args=([Histograms,Phases,restraintDict,rigidbodyDict],
                            parmDict,varyList,calcControls,pawleyLookup,dlg), 
                            nwalkers=nwalkers, nIter=nIterMCMC, Print=False)
        ndim=len(sig)     
        nToOmit = min(nIterMCMC/2, 50)                
        samples = sampler.chain[:, nToOmit:, :].reshape((-1, ndim))                                          
        bins = min(nIterMCMC*nwalkers/5, 20)
        fnametxt = ospath.splitext(GPXfile)[0]+'-MCMC.txt'
        np.savetxt(fnametxt, sampler.flatchain, header=str(varyList)) 
        try:
            import triangle
            fig = triangle.corner(samples, bins=bins, labels=varyList, 
                                  quantiles=[0.16, 0.5, 0.84])
            fnamepng = ospath.splitext(GPXfile)[0]+'-MCMC.png'
            fig.savefig(fnamepng)
        except:
            print ('Cannot allocate memory for the image. Image saved as', fnamepng)
        
        if 'Jacobian' in Controls['deriv type']:            
            result = so.leastsq(G2stMth.errRefine,values,Dfun=G2stMth.dervRefine,full_output=True,
                ftol=Ftol,col_deriv=True,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
        elif 'Hessian' in Controls['deriv type']:
            result = G2mth.HessianLSQ(G2stMth.errRefine,values,Hess=G2stMth.HessRefine,ftol=Ftol,maxcyc=maxCyc,Print=ifPrint,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
        else:           #'numeric'
            result = so.leastsq(G2stMth.errRefine,values,full_output=True,ftol=Ftol,epsfcn=1.e-8,factor=Factor,
                args=([Histograms,Phases,restraintDict,rigidbodyDict],parmDict,varyList,calcControls,pawleyLookup,dlg))
      						
    if IfOK:
        G2stMth.GetFobsSq(Histograms,Phases,parmDict,calcControls)
    return IfOK,Rvals,result,covMatrix,sig

def plotCorrections(nHist, histNames, E_mu, E_beta, nBlocks, x, cc_opt, bb_opt, dx_opt, ydiff, ystd, yexp, ycor, yuncor):
# MULTIPLICATIVE    
    if np.any(E_mu):
        fig=plt.figure()
        fig.canvas.set_window_title('Multiplicative Error Plots') 
        for i in range(nHist):
            plt.subplot(nHist, 1, i+1)
            lbl = histNames[i]+', E_mu='+str(E_mu[i])
            if i==0:
                plt.title('Multiplicative Error Estimation')
            plt.ylabel('Syst Error')  
            plt.xlabel(r'$\mathsf{2\theta}$')
            plt.plot(x[i], cc_opt[i], linestyle='--', label=lbl)
            plt.legend()
        f_manager = plt.get_current_fig_manager()
        if hasattr(f_manager.window, 'SetPosition'):
            f_manager.window.SetPosition((10,10))
        if hasattr(f_manager.window, 'move'):
            f_manager.window.move(10,10)   
        #plt.savefig("test1.png")
# ADDITIVE        
    if np.any(E_beta):
        fig=plt.figure()
        fig.canvas.set_window_title('Additive Error Plots') 
        for i in range(nHist):
            plt.subplot(nHist, 1, i+1)
            lbl = histNames[i]+', E_mu='+str(E_mu[i])
            if i==0:
                plt.title('Additive Error Estimation')
            plt.ylabel('Syst Error')  
            plt.xlabel(r'$\mathsf{2\theta}$')
            plt.plot(x[i], bb_opt[i], linestyle='--', label=lbl)
            plt.legend()
        f_manager = plt.get_current_fig_manager()
        if hasattr(f_manager.window, 'SetPosition'):
            f_manager.window.SetPosition((200,10))
        if hasattr(f_manager.window, 'move'):
            f_manager.window.move(200,10)     
        #plt.savefig("test2.png")
# PEAK_SHAPE        
    if np.any(nBlocks):
        fig=plt.figure()
        fig.canvas.set_window_title('Peak-shape Correction Plots') 
        for i in range(nHist):
            plt.subplot(nHist, 1, i+1)
            lbl = histNames[i]
            plt.xlabel(r'$\mathsf{2\theta}$')
            if i==0:
                plt.title('Peak-Shape Error Estimation')
            plt.ylabel(r'$\mathsf{\delta_x}$')  
            plt.plot(x[i], dx_opt[i], label=lbl)
            plt.legend()
        #plt.savefig("test3.png")                     
# RESIDUALS                        
    fig=plt.figure()
    fig.canvas.set_window_title('Residual Plots') 
    from scipy import stats
    for i in range(nHist):
        plt.subplot(nHist, 1, i+1)
        DD = stats.kstest(np.multiply(ydiff[i], ystd[i]), 'norm')
        lbl = histNames[i]+', D='+str(DD[0])
        if i==0:
            plt.title('Residuals')
        plt.ylabel('resid')  
        plt.xlabel(r'$\mathsf{2\theta}$')
        plt.plot(x[i], ydiff[i], linestyle='-', label=lbl)
        plt.plot(x[i], 2*ystd[i], linestyle='-', color="gray")
        plt.plot(x[i], -2*ystd[i], linestyle='-', color="gray")
        plt.legend()
    f_manager = plt.get_current_fig_manager()
    if hasattr(f_manager.window, 'SetPosition'):
        f_manager.window.SetPosition((400,400))
    if hasattr(f_manager.window, 'move'):
        f_manager.window.move(400,400)  
    #plt.savefig("test4.png")                   
# HISTOGRAMS
    fig=plt.figure()
    fig.canvas.set_window_title('HISTOGRAMS') 
    for i in range(nHist):
        plt.subplot(nHist, 1, i+1)
        lbl = histNames[i]
        if i==0:
            plt.title('HISTOGRAMS')
        plt.ylabel('Intensity')  
        plt.xlabel(r'$\mathsf{2\theta}$')
        plt.plot(x[i], yexp[i], 'ro', label=lbl)
        plt.plot(x[i], ycor[i], linestyle='-', color="red", label='optimized')
        plt.plot(x[i], yuncor[i], linestyle='-', color="gray", label='uncorrected')
        plt.legend()
    f_manager = plt.get_current_fig_manager()
    if hasattr(f_manager.window, 'SetPosition'):
        f_manager.window.SetPosition((800,200))
    if hasattr(f_manager.window, 'move'):
        f_manager.window.move(800,200)   
    #plt.savefig("test5.png", dpi = 600)
    plt.show()
# Anton Gagin />
     
	
def Refine(GPXfile,dlg=None,makeBack=True):
    'Global refinement -- refines to minimize against all histograms'
    import GSASIImpsubs as G2mp
    G2mp.InitMP()
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics

    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    G2stIO.ShowBanner(printFile)
    varyList = []
    parmDict = {}
    G2mv.InitVars()
    Controls = G2stIO.GetControls(GPXfile)
    G2stIO.ShowControls(Controls,printFile)
    calcControls = {}
    calcControls.update(Controls)
    constrDict,fixedList = G2stIO.GetConstraints(GPXfile)
    restraintDict = G2stIO.GetRestraints(GPXfile)
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print (' *** ERROR - you have no phases to refine! ***')
        print (' *** Refine aborted ***')
        return False,'No phases'
    if not Histograms:
        print (' *** ERROR - you have no data to refine with! ***')
        print (' *** Refine aborted ***')
        return False,'No data'
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,MFtables,maxSSwave = \
        G2stIO.GetPhaseData(Phases,restraintDict,rbIds,pFile=printFile)
    calcControls['atomIndx'] = atomIndx
    calcControls['Natoms'] = Natoms
    calcControls['FFtables'] = FFtables
    calcControls['BLtables'] = BLtables
    calcControls['MFtables'] = MFtables
    calcControls['maxSSwave'] = maxSSwave
    hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histograms,pFile=printFile)
    TwConstr,TwFixed = G2stIO.makeTwinFrConstr(Phases,Histograms,hapVary)
    constrDict += TwConstr
    fixedList += TwFixed
    calcControls.update(controlDict)
    histVary,histDict,controlDict = G2stIO.GetHistogramData(Histograms,pFile=printFile)
    calcControls.update(controlDict)
    varyList = rbVary+phaseVary+hapVary+histVary
    parmDict.update(rbDict)
    parmDict.update(phaseDict)
    parmDict.update(hapDict)
    parmDict.update(histDict)
    G2stIO.GetFprime(calcControls,Histograms)
    # do constraint processing
    varyListStart = tuple(varyList) # save the original varyList before dependent vars are removed
    try:
        groups,parmlist = G2mv.GroupConstraints(constrDict)
        G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList,parmDict)
        #print G2mv.VarRemapShow(varyList)
        #print 'DependentVars',G2mv.GetDependentVars()
        #print 'IndependentVars',G2mv.GetIndependentVars()
    except G2mv.ConstraintException:
        print (' *** ERROR - your constraints are internally inconsistent ***')
        #errmsg, warnmsg = G2mv.CheckConstraints(varyList,constrDict,fixedList)
        #print 'Errors',errmsg
        #if warnmsg: print 'Warnings',warnmsg
        return False,' Constraint error'
#    print G2mv.VarRemapShow(varyList)

    ifSeq = False
    printFile.write('\n Refinement results:\n')
    printFile.write(135*'-'+'\n')
    if True:
#    try:
        covData = {}
# </ Anton Gagin    		
        IfOK,Rvals,result,covMatrix,sig = RefineCore(Controls,Histograms,Phases,restraintDict,
            rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifSeq,printFile,dlg,GPXfile)
# Anton Gagin />			
        if IfOK:
            sigDict = dict(zip(varyList,sig))
            newCellDict = G2stMth.GetNewCellParms(parmDict,varyList)
            newAtomDict = G2stMth.ApplyXYZshifts(parmDict,varyList)
            covData = {'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                       'varyListStart':varyListStart,
                       'covMatrix':covMatrix,'title':GPXfile,'newAtomDict':newAtomDict,
                       'newCellDict':newCellDict,'freshCOV':True}
            # add the uncertainties into the esd dictionary (sigDict)
            sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))
            G2mv.PrintIndependentVars(parmDict,varyList,sigDict,pFile=printFile)
            G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
            G2stIO.SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
            G2stIO.SetPhaseData(parmDict,sigDict,Phases,rbIds,covData,restraintDict,printFile)
            G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histograms,calcControls['FFtables'],pFile=printFile)
            G2stIO.SetHistogramData(parmDict,sigDict,Histograms,calcControls['FFtables'],pFile=printFile)
            G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histograms,Phases,rigidbodyDict,covData,makeBack)
            printFile.close()
            print (' Refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst')
            print (' ***** Refinement successful *****')
        else:
            print ('****ERROR - Refinement failed')
            raise G2obj.G2Exception('****ERROR - Refinement failed')
#    except G2obj.G2Exception(Msg):
#        printFile.close()
#        return False,Msg.msg

#for testing purposes!!!
    if DEBUG and IfOK:
#needs: values,HistoPhases,parmDict,varylist,calcControls,pawleyLookup
        fl = open(ospath.splitext(GPXfile)[0]+'.testDeriv','wb')
        cPickle.dump(result[0],fl,1)
        cPickle.dump([Histograms,Phases,restraintDict,rigidbodyDict],fl,1)
        cPickle.dump([constrDict,fixedList,G2mv.GetDependentVars()],fl,1)
        cPickle.dump(parmDict,fl,1)
        cPickle.dump(varyList,fl,1)
        cPickle.dump(calcControls,fl,1)
        cPickle.dump(pawleyLookup,fl,1)
        fl.close()
    if dlg:
        return True,Rvals

def phaseCheck(phaseVary,Phases,histogram):
    '''
    Removes unused parameters from phase varylist if phase not in histogram
    #TODO - implement "Fix FXU" for seq refinement here - done?
    '''
    NewVary = []
    for phase in Phases:
        if histogram not in Phases[phase]['Histograms']: continue
        if Phases[phase]['Histograms'][histogram]['Use']:
            pId = Phases[phase]['pId']
            newVary = [item for item in phaseVary if item.split(':')[0] == str(pId)]
            FixVals = Phases[phase]['Histograms'][histogram].get('Fix FXU',' ')
            if 'F' in FixVals:
                newVary = [item for item in newVary if not 'Afrac' in item]
            if 'X' in FixVals:
                newVary = [item for item in newVary if not 'dA' in item]
            if 'U' in FixVals:
                newVary = [item for item in newVary if not 'AU' in item]
            if 'M' in FixVals:
                newVary = [item for item in newVary if not 'AM' in item]
            NewVary += newVary
    return NewVary

def SeqRefine(GPXfile,dlg,PlotFunction=None,G2frame=None):
    '''Perform a sequential refinement -- cycles through all selected histgrams,
    one at a time
    '''
    import GSASIImpsubs as G2mp
    G2mp.InitMP()
    import pytexture as ptx
    ptx.pyqlmninit()            #initialize fortran arrays for spherical harmonics

    printFile = open(ospath.splitext(GPXfile)[0]+'.lst','w')
    print ('Starting Sequential Refinement')
    G2stIO.ShowBanner(printFile)
    Controls = G2stIO.GetControls(GPXfile)
    G2stIO.ShowControls(Controls,printFile,SeqRef=True)
    restraintDict = G2stIO.GetRestraints(GPXfile)
    Histograms,Phases = G2stIO.GetUsedHistogramsAndPhases(GPXfile)
    if not Phases:
        print (' *** ERROR - you have no phases to refine! ***')
        print (' *** Refine aborted ***')
        return False,'No phases'
    if not Histograms:
        print (' *** ERROR - you have no data to refine with! ***')
        print (' *** Refine aborted ***')
        return False,'No data'
    rigidbodyDict = G2stIO.GetRigidBodies(GPXfile)
    rbIds = rigidbodyDict.get('RBIds',{'Vector':[],'Residue':[]})
    rbVary,rbDict = G2stIO.GetRigidBodyModels(rigidbodyDict,pFile=printFile)
    Natoms,atomIndx,phaseVary,phaseDict,pawleyLookup,FFtables,BLtables,MFtables,maxSSwave = \
        G2stIO.GetPhaseData(Phases,restraintDict,rbIds,False,printFile,seqRef=True)
    for item in phaseVary:
        if '::A0' in item:
            print ('**** WARNING - lattice parameters should not be refined in a sequential refinement ****')
            print ('****           instead use the Dij parameters for each powder histogram            ****')
            return False,'Lattice parameter refinement error - see console message'
        if '::C(' in item:
            print ('**** WARNING - phase texture parameters should not be refined in a sequential refinement ****')
            print ('****           instead use the C(L,N) parameters for each powder histogram               ****')
            return False,'Phase texture refinement error - see console message'
    if 'Seq Data' in Controls:
        histNames = Controls['Seq Data']
    else:
        histNames = G2stIO.GetHistogramNames(GPXfile,['PWDR',])
    if Controls.get('Reverse Seq'):
        histNames.reverse()
    SeqResult = G2stIO.GetSeqResult(GPXfile)
#    SeqResult = {'SeqPseudoVars':{},'SeqParFitEqList':[]}
    makeBack = True
    Histo = {}
    NewparmDict = {}
    for ihst,histogram in enumerate(histNames):
        print('\nRefining with '+str(histogram))
        ifPrint = False
        if dlg:
            dlg.SetTitle('Residual for histogram '+str(ihst))
        calcControls = {}
        calcControls['atomIndx'] = atomIndx
        calcControls['Natoms'] = Natoms
        calcControls['FFtables'] = FFtables
        calcControls['BLtables'] = BLtables
        calcControls['MFtables'] = MFtables
        calcControls['maxSSwave'] = maxSSwave
        if histogram not in Histograms:
            print("Error: not found!")
            continue
    #TODO - implement "Fix FXU" for seq refinement here - done?
        hId = Histograms[histogram]['hId']
        redphaseVary = phaseCheck(phaseVary,Phases,histogram)
        Histo = {histogram:Histograms[histogram],}
        hapVary,hapDict,controlDict = G2stIO.GetHistogramPhaseData(Phases,Histo,Print=False)
        calcControls.update(controlDict)
        histVary,histDict,controlDict = G2stIO.GetHistogramData(Histo,False)
        calcControls.update(controlDict)
        varyList = rbVary+redphaseVary+hapVary+histVary
#        if not ihst:
            # save the initial vary list, but without histogram numbers on parameters
        saveVaryList = varyList[:]
        for i,item in enumerate(saveVaryList):
            items = item.split(':')
            if items[1]:
                items[1] = ''
            item = ':'.join(items)
            saveVaryList[i] = item
        if not ihst:
            SeqResult['varyList'] = saveVaryList
        else:
            SeqResult['varyList'] = list(set(SeqResult['varyList']+saveVaryList))
        parmDict = {}
        parmDict.update(phaseDict)
        parmDict.update(hapDict)
        parmDict.update(histDict)
        if Controls['Copy2Next']:
            #parmDict.update(NewparmDict) # don't use in case extra entries would cause a problem
            for parm in NewparmDict:
                if parm in parmDict:
                    parmDict[parm] = NewparmDict[parm]
        G2stIO.GetFprime(calcControls,Histo)
        # do constraint processing
        #reload(G2mv) # debug
        G2mv.InitVars()
        constrDict,fixedList = G2stIO.GetConstraints(GPXfile)
        varyListStart = tuple(varyList) # save the original varyList before dependent vars are removed
        try:
            groups,parmlist = G2mv.GroupConstraints(constrDict)
            G2mv.GenerateConstraints(groups,parmlist,varyList,constrDict,fixedList,parmDict,SeqHist=hId)
#            if GSASIIpath.GetConfigValue('debug'): print("DBG_"+
#                G2mv.VarRemapShow(varyList,True))
            constraintInfo = (groups,parmlist,constrDict,fixedList,ihst)
        except G2mv.ConstraintException:
            print (' *** ERROR - your constraints are internally inconsistent ***')
            #errmsg, warnmsg = G2mv.CheckConstraints(varyList,constrDict,fixedList)
            #print 'Errors',errmsg
            #if warnmsg: print 'Warnings',warnmsg
            return False,' Constraint error'
        #print G2mv.VarRemapShow(varyList)
        if not ihst:
            # first histogram to refine against
            firstVaryList = []
            for item in varyList:
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                firstVaryList.append(item)
            newVaryList = firstVaryList
        else:
            newVaryList = []
            for item in varyList:
                items = item.split(':')
                if items[1]:
                    items[1] = ''
                item = ':'.join(items)
                newVaryList.append(item)
        if newVaryList != firstVaryList and Controls['Copy2Next']:
            # variable lists are expected to match between sequential refinements when Copy2Next is on
            #print '**** ERROR - variable list for this histogram does not match previous'
            #print '     Copy of variables is not possible'
            #print '\ncurrent histogram',histogram,'has',len(newVaryList),'variables'
            combined = list(set(firstVaryList+newVaryList))
            c = [var for var in combined if var not in newVaryList]
            p = [var for var in combined if var not in firstVaryList]
            print('*** Variables change ***')
            for typ,vars in [('Removed',c),('Added',p)]:
                line = '  '+typ+': '
                if vars:
                    for var in vars:
                        if len(line) > 70:
                            print(line)
                            line = '    '
                        line += var + ', '
                else:
                        line += 'none, '
                print(line[:-2])
            firstVaryList = newVaryList

        ifSeq = True
        printFile.write('\n Refinement results for histogram: %s\n'%histogram)
        printFile.write(135*'-'+'\n')
        if True:
#        try:
            IfOK,Rvals,result,covMatrix,sig = RefineCore(Controls,Histo,Phases,restraintDict,
                rigidbodyDict,parmDict,varyList,calcControls,pawleyLookup,ifSeq,printFile,dlg)
            if PlotFunction:
                PlotFunction(G2frame,Histo[histogram]['Data'],histogram)

            print ('  wR = %7.2f%%, chi**2 = %12.6g, reduced chi**2 = %6.2f, last delta chi = %.4f'%(
                Rvals['Rwp'],Rvals['chisq'],Rvals['GOF']**2,Rvals['DelChi2']))
            # add the uncertainties into the esd dictionary (sigDict)
            if not IfOK:
                print('***** Sequential refinement failed at histogram '+histogram)
                break
            sigDict = dict(zip(varyList,sig))
            # the uncertainties for dependent constrained parms into the esd dict
            sigDict.update(G2mv.ComputeDepESD(covMatrix,varyList,parmDict))

            # a dict with values & esds for dependent (constrained) parameters - avoid extraneous holds
            depParmDict = {i:(parmDict[i],sigDict[i]) for i in varyListStart if i in sigDict and i not in varyList}
            newCellDict = copy.deepcopy(G2stMth.GetNewCellParms(parmDict,varyList))
            newAtomDict = copy.deepcopy(G2stMth.ApplyXYZshifts(parmDict,varyList))
            histRefData = {
                'variables':result[0],'varyList':varyList,'sig':sig,'Rvals':Rvals,
                'varyListStart':varyListStart,
                'covMatrix':covMatrix,'title':histogram,'newAtomDict':newAtomDict,
                'newCellDict':newCellDict,'depParmDict':depParmDict,
                'constraintInfo':constraintInfo,
                'parmDict':parmDict}
            SeqResult[histogram] = histRefData
            G2stMth.ApplyRBModels(parmDict,Phases,rigidbodyDict,True)
    #        G2stIO.SetRigidBodyModels(parmDict,sigDict,rigidbodyDict,printFile)
            G2stIO.SetHistogramPhaseData(parmDict,sigDict,Phases,Histo,None,ifPrint,printFile)
            G2stIO.SetHistogramData(parmDict,sigDict,Histo,None,ifPrint,printFile)
            G2stIO.SetUsedHistogramsAndPhases(GPXfile,Histo,Phases,rigidbodyDict,histRefData,makeBack)
            makeBack = False
            NewparmDict = {}
            # make dict of varied parameters in current histogram, renamed to
            # next histogram, for use in next refinement.
            if Controls['Copy2Next'] and ihst < len(histNames)-1:
                hId = Histo[histogram]['hId'] # current histogram
                nexthId = Histograms[histNames[ihst+1]]['hId']
                for parm in set(list(varyList)+list(varyListStart)):
                    items = parm.split(':')
                    if len(items) < 3: continue
                    if str(hId) in items[1]:
                        items[1] = str(nexthId)
                        newparm = ':'.join(items)
                        NewparmDict[newparm] = parmDict[parm]
#        except G2obj.G2Exception(Msg):
#            printFile.close()
#            print (' ***** Refinement aborted *****')
#            return False,Msg.msg
    SeqResult['histNames'] = [itm for itm in G2stIO.GetHistogramNames(GPXfile,['PWDR',]) if itm in SeqResult.keys()]
    G2stIO.SetSeqResult(GPXfile,Histograms,SeqResult)
    printFile.close()
    print (' Sequential refinement results are in file: '+ospath.splitext(GPXfile)[0]+'.lst')
    print (' ***** Sequential refinement successful *****')
    return True,'Success'

def RetDistAngle(DisAglCtls,DisAglData):
    '''Compute and return distances and angles

    :param dict DisAglCtls: contains distance/angle radii usually defined using
       :func:`GSASIIctrlGUI.DisAglDialog`
    :param dict DisAglData: contains phase data:
       Items 'OrigAtoms' and 'TargAtoms' contain the atoms to be used
       for distance/angle origins and atoms to be used as targets.
       Item 'SGData' has the space group information (see :ref:`Space Group object<SGData_table>`)

    :returns: AtomLabels,DistArray,AngArray where:

      **AtomLabels** is a dict of atom labels, keys are the atom number

      **DistArray** is a dict keyed by the origin atom number where the value is a list
      of distance entries. The value for each distance is a list containing:

        0) the target atom number (int);
        1) the unit cell offsets added to x,y & z (tuple of int values)
        2) the symmetry operator number (which may be modified to indicate centering and center of symmetry)
        3) an interatomic distance in A (float)
        4) an uncertainty on the distance in A or 0.0 (float)

      **AngArray** is a dict keyed by the origin (central) atom number where
      the value is a list of
      angle entries. The value for each angle entry consists of three values:

        0) a distance item reference for one neighbor (int)
        1) a distance item reference for a second neighbor (int)
        2) a angle, uncertainty pair; the s.u. may be zero (tuple of two floats)

      The AngArray distance reference items refer directly to the index of the items in the
      DistArray item for the list of distances for the central atom.
    '''
    import numpy.ma as ma

    SGData = DisAglData['SGData']
    Cell = DisAglData['Cell']

    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    if 'covData' in DisAglData:
        covData = DisAglData['covData']
        covMatrix = covData['covMatrix']
        varyList = covData['varyList']
        pfx = str(DisAglData['pId'])+'::'

    Factor = DisAglCtls['Factors']
    Radii = dict(zip(DisAglCtls['AtomTypes'],zip(DisAglCtls['BondRadii'],DisAglCtls['AngleRadii'])))
    indices = (-2,-1,0,1,2)
    Units = np.array([[h,k,l] for h in indices for k in indices for l in indices])
    origAtoms = DisAglData['OrigAtoms']
    targAtoms = DisAglData['TargAtoms']
    AtomLabels = {}
    for Oatom in origAtoms:
        AtomLabels[Oatom[0]] = Oatom[1]
    for Oatom in targAtoms:
        AtomLabels[Oatom[0]] = Oatom[1]
    DistArray = {}
    AngArray = {}
    for Oatom in origAtoms:
        DistArray[Oatom[0]] = []
        AngArray[Oatom[0]] = []
        OxyzNames = ''
        IndBlist = []
        Dist = []
        Vect = []
        VectA = []
        angles = []
        for Tatom in targAtoms:
            Xvcov = []
            TxyzNames = ''
            if 'covData' in DisAglData:
                OxyzNames = [pfx+'dAx:%d'%(Oatom[0]),pfx+'dAy:%d'%(Oatom[0]),pfx+'dAz:%d'%(Oatom[0])]
                TxyzNames = [pfx+'dAx:%d'%(Tatom[0]),pfx+'dAy:%d'%(Tatom[0]),pfx+'dAz:%d'%(Tatom[0])]
                Xvcov = G2mth.getVCov(OxyzNames+TxyzNames,varyList,covMatrix)
            result = G2spc.GenAtom(Tatom[3:6],SGData,False,Move=False)
            BsumR = (Radii[Oatom[2]][0]+Radii[Tatom[2]][0])*Factor[0]
            AsumR = (Radii[Oatom[2]][1]+Radii[Tatom[2]][1])*Factor[1]
            for [Txyz,Top,Tunit,Spn] in result:
                Dx = (Txyz-np.array(Oatom[3:6]))+Units
                dx = np.inner(Amat,Dx)
                dist = ma.masked_less(np.sqrt(np.sum(dx**2,axis=0)),0.5)
                IndB = ma.nonzero(ma.masked_greater(dist-BsumR,0.))
                if np.any(IndB):
                    for indb in IndB:
                        for i in range(len(indb)):
                            if str(dx.T[indb][i]) not in IndBlist:
                                IndBlist.append(str(dx.T[indb][i]))
                                unit = Units[indb][i]
                                tunit = (unit[0]+Tunit[0],unit[1]+Tunit[1],unit[2]+Tunit[2])
                                pdpx = G2mth.getDistDerv(Oatom[3:6],Tatom[3:6],Amat,unit,Top,SGData)
                                sig = 0.0
                                if len(Xvcov):
                                    sig = np.sqrt(np.inner(pdpx,np.inner(pdpx,Xvcov)))
                                Dist.append([Oatom[0],Tatom[0],tunit,Top,ma.getdata(dist[indb])[i],sig])
                                if (Dist[-1][-2]-AsumR) <= 0.:
                                    Vect.append(dx.T[indb][i]/Dist[-1][-2])
                                    VectA.append([OxyzNames,np.array(Oatom[3:6]),TxyzNames,np.array(Tatom[3:6]),unit,Top])
                                else:
                                    Vect.append([0.,0.,0.])
                                    VectA.append([])
        for D in Dist:
            DistArray[Oatom[0]].append(D[1:])
        Vect = np.array(Vect)
        angles = np.zeros((len(Vect),len(Vect)))
        angsig = np.zeros((len(Vect),len(Vect)))
        for i,veca in enumerate(Vect):
            if np.any(veca):
                for j,vecb in enumerate(Vect):
                    if np.any(vecb):
                        angles[i][j],angsig[i][j] = G2mth.getAngSig(VectA[i],VectA[j],Amat,SGData,covData)
                        if i <= j: continue
                        AngArray[Oatom[0]].append((i,j,
                            G2mth.getAngSig(VectA[i],VectA[j],Amat,SGData,covData)))
    return AtomLabels,DistArray,AngArray

def PrintDistAngle(DisAglCtls,DisAglData,out=sys.stdout):
    '''Print distances and angles

    :param dict DisAglCtls: contains distance/angle radii usually defined using
       :func:`GSASIIctrlGUI.DisAglDialog`
    :param dict DisAglData: contains phase data:
       Items 'OrigAtoms' and 'TargAtoms' contain the atoms to be used
       for distance/angle origins and atoms to be used as targets.
       Item 'SGData' has the space group information (see :ref:`Space Group object<SGData_table>`)
    :param file out: file object for output. Defaults to sys.stdout.
    '''
    def MyPrint(s):
        out.write(s+'\n')
        # print(s,file=out) # use in Python 3

    def ShowBanner(name):
        MyPrint(80*'*')
        MyPrint('   Interatomic Distances and Angles for phase '+name)
        MyPrint((80*'*')+'\n')

    ShowBanner(DisAglCtls['Name'])
    SGData = DisAglData['SGData']
    SGtext,SGtable = G2spc.SGPrint(SGData)
    for line in SGtext: MyPrint(line)
    if len(SGtable) > 1:
        for i,item in enumerate(SGtable[::2]):
            if 2*i+1 == len(SGtable):
                line = ' %s'%(item.ljust(30))
            else:
                line = ' %s %s'%(item.ljust(30),SGtable[2*i+1].ljust(30))
            MyPrint(line)
    else:
        MyPrint(' ( 1)    %s'%(SGtable[0])) #triclinic case
    Cell = DisAglData['Cell']

    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    if 'covData' in DisAglData:
        covData = DisAglData['covData']
        pfx = str(DisAglData['pId'])+'::'
        A = G2lat.cell2A(Cell[:6])
        cellSig = G2stIO.getCellEsd(pfx,SGData,A,covData)
        names = [' a = ',' b = ',' c = ',' alpha = ',' beta = ',' gamma = ',' Volume = ']
        valEsd = [G2mth.ValEsd(Cell[i],cellSig[i],True) for i in range(7)]
        line = '\n Unit cell:'
        for name,vals in zip(names,valEsd):
            line += name+vals
        MyPrint(line)
    else:
        MyPrint('\n Unit cell: a = '+('%.5f'%Cell[0])+' b = '+('%.5f'%Cell[1])+' c = '+('%.5f'%Cell[2])+
            ' alpha = '+('%.3f'%Cell[3])+' beta = '+('%.3f'%Cell[4])+' gamma = '+
            ('%.3f'%Cell[5])+' Volume = '+('%.3f'%Cell[6]))

    AtomLabels,DistArray,AngArray = RetDistAngle(DisAglCtls,DisAglData)
    origAtoms = DisAglData['OrigAtoms']
    for Oatom in origAtoms:
        i = Oatom[0]
        Dist = DistArray[i]
        nDist = len(Dist)
        angles = np.zeros((nDist,nDist))
        angsig = np.zeros((nDist,nDist))
        for k,j,tup in AngArray[i]:
            angles[k][j],angsig[k][j] = angles[j][k],angsig[j][k] = tup
        line = ''
        for i,x in enumerate(Oatom[3:6]):
            line += ('%12.5f'%x).rstrip('0')
        MyPrint('\n Distances & angles for '+Oatom[1]+' at '+line.rstrip())
        MyPrint(80*'*')
        line = ''
        for dist in Dist[:-1]:
            line += '%12s'%(AtomLabels[dist[0]].center(12))
        MyPrint('  To       cell +(sym. op.)      dist.  '+line.rstrip())
        for i,dist in enumerate(Dist):
            line = ''
            for j,angle in enumerate(angles[i][0:i]):
                sig = angsig[i][j]
                if angle:
                    if sig:
                        line += '%12s'%(G2mth.ValEsd(angle,sig,True).center(12))
                    else:
                        val = '%.3f'%(angle)
                        line += '%12s'%(val.center(12))
                else:
                    line += 12*' '
            if dist[4]:            #sig exists!
                val = G2mth.ValEsd(dist[3],dist[4])
            else:
                val = '%8.4f'%(dist[3])
            tunit = '[%2d%2d%2d]'% dist[1]
            MyPrint(('  %8s%10s+(%4d) %12s'%(AtomLabels[dist[0]].ljust(8),tunit.ljust(10),dist[2],val.center(12)))+line.rstrip())

def DisAglTor(DATData):
    'Needs a doc string'
    SGData = DATData['SGData']
    Cell = DATData['Cell']

    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    covData = {}
    pfx = ''
    if 'covData' in DATData:
        covData = DATData['covData']
        pfx = str(DATData['pId'])+'::'
    Datoms = []
    Oatoms = []
    for i,atom in enumerate(DATData['Datoms']):
        symop = atom[-1].split('+')
        if len(symop) == 1:
            symop.append('0,0,0')
        symop[0] = int(symop[0])
        symop[1] = eval(symop[1])
        atom.append(symop)
        Datoms.append(atom)
        oatom = DATData['Oatoms'][i]
        names = ['','','']
        if pfx:
            names = [pfx+'dAx:'+str(oatom[0]),pfx+'dAy:'+str(oatom[0]),pfx+'dAz:'+str(oatom[0])]
        oatom += [names,]
        Oatoms.append(oatom)
    atmSeq = [atom[1]+'('+atom[-2]+')' for atom in Datoms]
    if DATData['Natoms'] == 4:  #torsion
        Tors,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        print (' Torsion angle for %s atom sequence: %s = %s'%(DATData['Name'],str(atmSeq).replace("'","")[1:-1],G2mth.ValEsd(Tors,sig)))
        print (' NB: Atom sequence determined by selection order')
        return      # done with torsion
    elif DATData['Natoms'] == 3:  #angle
        Ang,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        print (' Angle in %s for atom sequence: %s = %s'%(DATData['Name'],str(atmSeq).replace("'","")[1:-1],G2mth.ValEsd(Ang,sig)))
        print (' NB: Atom sequence determined by selection order')
    else:   #2 atoms - distance
        Dist,sig = G2mth.GetDATSig(Oatoms,Datoms,Amat,SGData,covData)
        print (' Distance in %s for atom sequence: %s = %s'%(DATData['Name'],str(atmSeq).replace("'","")[1:-1],G2mth.ValEsd(Dist,sig)))

def BestPlane(PlaneData):
    'Needs a doc string'

    def ShowBanner(name):
        print (80*'*')
        print ('   Best plane result for phase '+name)
        print (80*'*','\n')

    ShowBanner(PlaneData['Name'])

    Cell = PlaneData['Cell']
    Amat,Bmat = G2lat.cell2AB(Cell[:6])
    Atoms = PlaneData['Atoms']
    sumXYZ = np.zeros(3)
    XYZ = []
    Natoms = len(Atoms)
    for atom in Atoms:
        xyz = np.array(atom[3:6])
        XYZ.append(xyz)
        sumXYZ += xyz
    sumXYZ /= Natoms
    XYZ = np.array(XYZ)-sumXYZ
    XYZ = np.inner(Amat,XYZ).T
    Zmat = np.zeros((3,3))
    for i,xyz in enumerate(XYZ):
        Zmat += np.outer(xyz.T,xyz)
    print (' Selected atoms centered at %10.5f %10.5f %10.5f'%(sumXYZ[0],sumXYZ[1],sumXYZ[2]))
    Evec,Emat = nl.eig(Zmat)
    Evec = np.sqrt(Evec)/(Natoms-3)
    Order = np.argsort(Evec)
    XYZ = np.inner(XYZ,Emat.T).T
    XYZ = np.array([XYZ[Order[2]],XYZ[Order[1]],XYZ[Order[0]]]).T
    print (' Atoms in Cartesian best plane coordinates:')
    print (' Name         X         Y         Z')
    for i,xyz in enumerate(XYZ):
        print (' %6s%10.3f%10.3f%10.3f'%(Atoms[i][1].ljust(6),xyz[0],xyz[1],xyz[2]))
    print ('\n Best plane RMS X =%8.3f, Y =%8.3f, Z =%8.3f'%(Evec[Order[2]],Evec[Order[1]],Evec[Order[0]]))

def main():
    'Called to run a refinement when this module is executed '
    starttime = time.time()
    arg = sys.argv
    if len(arg) > 1:
        GPXfile = arg[1]
        if not ospath.exists(GPXfile):
            print ('ERROR - '+GPXfile+" doesn't exist!")
            exit()
    else:
        print ('ERROR - missing filename')
        exit()
    # TODO: figure out if this is a sequential refinement and call SeqRefine(GPXfile,None)
    Refine(GPXfile,None)
    print("Done. Execution time {:.2f} sec.".format(time.time()-starttime))

if __name__ == '__main__':
    GSASIIpath.InvokeDebugOpts()
    main()
