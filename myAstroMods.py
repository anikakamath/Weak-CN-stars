# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 12:06:04 2016

@author: jonhays
"""


import numpy as np
import matplotlib.pyplot as py

#%%
def CMD_sample_stars(data, NUM_SAMPLES, s2nThresh=10, magThresh=19):
          
    b2ir = data['F475W'] - data['F814W']                                        # Array of blue to infrared ratio (temp) of non-carbon stars
    
    maxb2ir = np.nanmax(b2ir)                                                   # Finding max value of B - I (blue to I.R. ratio)
    minb2ir = np.nanmin(b2ir)                                                   # Finding min value of B - I

    upperLimit = 14.63700066                                                    # This upper limit was determined to be the highest value for B - I for reasonable data. Anythinging above this value was due to values for F475W being set to 99.99900055 
    binLimits = np.linspace(minb2ir, upperLimit, NUM_SAMPLES + 1)               # by default when no reading occurs.       
    
    binIndices = []                         
    for i in range(0, NUM_SAMPLES):                                         
        tempArray = []                          
        for j in range(0, len(b2ir)):                
            if b2ir[j] >= binLimits[i] and b2ir[j] <= binLimits[i + 1]:
                tempArray += [j]
                    
        binIndices += [tempArray] 
        
                                                                                
    maxS2N = np.zeros(len(binIndices))                                         # Finds highest signal to noise for individual 
    sampleIndices = np.zeros(len(binIndices))                                  # and place corresponding indices in an array with NUM_SAMPLES elements.
    for i in range(0, len(binIndices)):                                        # (if non-empty)
        if binIndices[i]:
            maxS2N[i] = data['ZSNR'][binIndices[i][0]]
            sampleIndices[i] = binIndices[i][0]
            for j in range(1, len(binIndices[i])):
                if maxS2N[i] < data['ZSNR'][binIndices[i][j]]:
                    maxS2N[i] = data['ZSNR'][binIndices[i][j]]
                    sampleIndices[i] = binIndices[i][j]
                else:
                    pass
        else:
            maxS2N[i] = np.nan
            sampleIndices[i] = np.nan
            
    initIndices = sampleIndices   
    sampleIndices = np.array(sampleIndices[np.isfinite(sampleIndices)])
    
    #    sampleIndices = [i for i in sampleIndices if np.isfinite(data['F475W'][i]) and np.isfinite(data['F814W'][i])]                  # removes NaN values from ARRAY (must be array to use this built in numpy function)
    sampleIndices = [i for i in sampleIndices if data['ZSNR'][i] > s2nThresh and data['F814W'][i] < magThresh and data['F814w'][i] > 16 and data['ZQUAL'][i] != 2 and data['MASK'][i][0:4] == 'mct6' and data['MASK'][i][4] not in ['A', 'J', 'N', 'Y', 'Z']]  # This line removes any indices with corresponding S/N values that are less than s2nThresh (default value is 10),.
    sampleData = data[sampleIndices]                                                                               # as well as indices with corresponding magnitude values that are higher than magThresh (default value is 15).


    
    indicesRemoved = len(initIndices) - len(sampleIndices)
    indicesRemaining = len(sampleIndices)
    
    print('Returned a total of %d samples; %d samples were removed due to low S/N values, high magnitude values (low brightness), NaN values, and stars with a ZQUAL value of 2.' % (indicesRemaining, indicesRemoved) )    
    
    return [sampleData, sampleIndices]

  

#%%

def weighted_sum_array(scData, strongCindex, sampleData):
    
    scSpec = scData['SPEC'][strongCindex]                                      # Spectral data for the strong carbon 'test' star
    wArray = [[.1, .9], [.2, .8], [.3, .7], [.4, .6], [.5, .5], [.6, .4], [.7, .3], [.8, .2], [.9, .1]]  # array of possible wieghting combos
    specSize = 9230
    allWeightedSpectra = []
    for r in range(0, len(sampleData)):  
        sampleSpec = sampleData['SPEC'][r]                                                              # number of flux values in each spectrum array
        weightedSpectra = []
        tempMatrix = []
        for i in range(0, len(wArray)):
            tempArray = []
            for j in range(0, specSize):
                tempArray += [(wArray[i][0]*scSpec[j] + wArray[i][1]*sampleSpec[j])]
            weightedSpectra += [tempArray]
            tempMatrix += [weightedSpectra]
        allWeightedSpectra += [tempMatrix]
    
    return allWeightedSpectra
    
#%%
def weighted_sum(spec1, spec2, w1, w2):
    weighted_sum = w1*spec1 + w2*spec2
    return weighted_sum

#%%

def normal(specData, varData=None, lowerThresh=6924, upperThresh=7695):                       # 6924 and 7695 are indices of the spectral array that correspond to 8500 and 9000 respectively..
    centralSpec = [specData[i] for i in range(0, len(specData)) if i <= upperThresh and i >= lowerThresh]  
    centralSpec = centralSpec[centralSpec != np.nan]
    nor = 1 / np.median(centralSpec)
    normalSpec = nor*specData
    ivarNormalize = (1/nor)**2
    if varData is not None:
        norVar = ivarNormalize*varData
        return [normalSpec, norVar]
    else:
        return normalSpec

#%%
def normal_set(data, upperThresh=9000, lowerThresh=8500):  
  
    normalSpecSet = []
    nor = []
    for i in range(0, len(data)): 
        centralSpec = []          
        for j in range(0, len(data['SPEC'][i])):
            if (data['LBIN'][i][j] <= upperThresh) and (data['LBIN'][i][j] >= lowerThresh):
                centralSpec += [data['SPEC'][i][j]]
        nor += [1/np.nanmedian(centralSpec)]
        normalSpecSet += [nor[i]*data['SPEC'][i]]
        
    return [normalSpecSet, nor]

#%%

def cmd_plot(data, indices, s=20, clr='b', marker='o'):

    py.rcParams['figure.figsize'] = 11, 10
    py.scatter((data.F475W[indices] - data.F814W[indices]), data.F814W[indices], s=s, color=clr, marker=marker)
    py.xlim(-2, 8)
    py.ylim(24, 16)
    py.plot([1, 3], [24, 16], color='k', linestyle='--', linewidth=2)
    py.plot([-2, 10], [21.5, 21.5], color='k', linestyle='--', linewidth=2)  
    
def cmd_plot2(data, indices, s=20, clr='b', marker='o', label = None):
    py.scatter((data.F814W[indices] - data.F160W[indices]), data.F160W[indices], s=s, color=clr, marker=marker, label = label)

  

 #%%

def plotspectra(data, indices, ymax=25000, ysize=25):
        
    py.rcParams['figure.figsize'] = 20, ysize
    for i in range(0, len(indices)):
        if i % 3 == 0:
            clr = 'b'
        elif i % 3 == 1: clr = 'g'
        else: clr = 'r'
        py.plot(data['LBIN'][indices[i]], data['SPEC'][indices[i]] + 4000*i, color = clr)
    py.xlim(5000, 10000)
    py.ylim(-400, ymax)
    
def graph_legend(param):
    py.plot([1,2,3], label="Non-carbon stars (reference sample)", color = 'w')
    py.plot([3,2,1], label="Extreme carbon stars", color= 'w')
    py.plot([1,2,3], label="Strong carbon stars",color= 'w')
    py.plot([1,2,3], label="Medium carbon stars", color= 'w')
    py.plot([1,2,3], label="Weak carbon stars", color= 'w')
    py.plot([1,2,3], label="Marginal carbon stars", color= 'w')
    py.plot([1,2,3], label="Weak/marginal carbon bin", color= 'w')
    py.plot([1,2,3], label="Strong carbon bin", color= 'w')
        # Place a legend to the right of this smaller figure.
    py.legend(bbox_to_anchor=(param, 1),loc=2,fontsize=10,fancybox=True).get_frame().set_alpha(0.5)


    py.show()
    
#%%
def sigclip_coadd_arr(lbin, spec, ivar, **kwargs):
    '''
    Use this version if you want to input arrays of wavelengh, flux, and inverse variance
    Keywords are:
       sigLim = at what sigma level do you want your coadd clipped. Default is 3.5
       verbose = if True, output coadded as well as coadded+clipped spectra. Default is False.
    '''
    
    sigLim = kwargs.get('sigLim',3.5)
    verbose = kwargs.get('verbose',False)

    nlbin = len(lbin[0])
    nstars = len(lbin)

    coaddivar = np.nansum(ivar, axis = 0)
    coaddspec = np.nansum(spec*ivar, axis = 0)/coaddivar

    clipivar = []
    subspec = []

    for j in range(nlbin):
        if coaddivar[j] == 0:
            clipivar.append(0)
            subspec.append(np.nan)
        else:
            s = spec[:,j]
            iv = ivar[:,j]
            nsig = np.abs((s - coaddspec[j])*np.sqrt(iv))
            subspec.append(np.nansum(s[nsig < sigLim]*iv[nsig < sigLim]))
            clipivar.append(np.nansum(iv[nsig < sigLim]))

    clipspec = np.array(subspec)/np.array(clipivar)
    clipivar = np.array(clipivar)

    if verbose:
        return coaddspec, coaddivar, clipspec, clipivar
    else:
        return clipspec, clipivar

def sigclip_coadd_med(lbin, spec, ivar, **kwargs):
    '''
    MODIFIED: This version uses median instead of weighted sum to calculate degree of deviation for individual spectra.
    Use this version if you want to input arrays of wavelengh, flux, and inverse variance
    Keywords are:
       sigLim = at what sigma level do you want your coadd clipped. Default is 3.5
       verbose = if True, output coadded as well as coadded+clipped spectra. Default is False.
    '''
    
    sigLim = kwargs.get('sigLim',3.5)
    verbose = kwargs.get('verbose',False)

    nlbin = len(lbin[0])                                 # number of pixels
    nstars = len(lbin)                                   # number of stars

    coaddivar = np.nansum(ivar, axis = 0)                # summing up corresponding elements of ALL ivar arrays
    specMedian = np.nanmedian(spec, axis = 0) 
    clipivar = []                                        # Initialize empty list
    subspec = []

    for j in range(nlbin):                                # For each pixel,
        if coaddivar[j] == 0:                             # if the ivar value for said pixel == 0 for every ivar array,
            clipivar.append(0)                            # then clipivar[j] = 0, and subspec[j] = np.nan.               
            subspec.append(np.nan)                        # (.append adds value in parentheses to end of existing list)
        else:
            s = spec[:,j]                                 # s is an array that contains flux values for every star at pixel j 
            iv = ivar[:,j]                                # iv is an array that contains ivar values for every star at pixel j
            nsig = np.abs((s - specMedian[j])*np.sqrt(iv)) # subtract summed flux value at pixel j from each flux value in s, 
            subspec.append(np.nansum(s[nsig < sigLim]*iv[nsig < sigLim])) # then divide by standard deviation (since... 
            clipivar.append(np.nansum(iv[nsig < sigLim]))

    clipspec = np.array(subspec)/np.array(clipivar)
    clipivar = np.array(clipivar)

    if verbose:
        return specMedian, coaddivar, clipspec, clipivar
    else:
        return clipspec, clipivar
    
#%%

def coadd_med(fluxnorm, ivarnorm, lamb, siglim=3., verbose=False):
    
    '''
    Fluxnorm and ivarnorm are arrays of normalized flux and ivar values for multiple stars that will be coadded.
    lamb is a single array of wavelength values that corresponds to all the stars, so that len(lamb) == # of pixels.
    This assumes that every star in the coadd has the same wavelength array corresponding to its flux values.
    '''
    siglim = float(siglim)
    medflux=np.nanmedian(fluxnorm,axis=0)
    sigdiff=abs(fluxnorm-medflux)*np.sqrt(ivarnorm)
    clip=np.ones((len(fluxnorm), len(lamb)))
    clip[sigdiff>siglim]=0

    coadd=np.nansum((fluxnorm*ivarnorm)*clip,0)/np.nansum(ivarnorm*clip,0)
    coadderr=1/np.sqrt(np.nansum(ivarnorm*clip,0))
    
    if verbose:
        return coadd, coadderr
    else:
        return coadd
    

def calcSN_norm(l,f,**kwargs):
    #Calculates the S/N of a normalized (FLAT!) spectrum
    #Required Arguments:
    #  - l: wavelength (units irrelevant)
    #  - f: flux (units irrelevant)
    #  - ivar = ivar OR var = var
      
    v = kwargs.get('var',None)
    i = kwargs.get('ivar',None)

    if (v is None) & (i is None):
        print ('Needs ivar or variance')
        return np.nan
    
    if (v is not None) & (i is None):
        i = 1./v

    fGood = f[i > 0]
    iGood = i[i > 0]

    if len(fGood) > 1:
        sn = np.median(fGood)*np.sqrt(np.median(iGood))
    else:
        print ('Too few good pixels')
        sn = np.nan

    return sn
    
# def fit_spec(scienceSpec, testSpectra, lamb, **kwargs)