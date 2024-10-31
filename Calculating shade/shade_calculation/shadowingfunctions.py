# -*- coding: utf-8 -*-
# Ready for python action!
import numpy as np
# import matplotlib.pylab as plt
# from numba import jit

def shadowingfunctionglobalradiation(a, azimuth, altitude, scale, forsvf):

    #%This m.file calculates shadows on a DEM
    #% conversion
    degrees = np.pi/180.
    # if azimuth == 0.0:
        # azimuth = 0.000000000001
    azimuth = np.dot(azimuth, degrees)
    altitude = np.dot(altitude, degrees)
    #% measure the size of the image
    sizex = a.shape[0]
    sizey = a.shape[1]
    # if forsvf == 0:
    #     barstep = np.max([sizex, sizey])
    #     # dlg.progressBar.setRange(0, barstep)
    #% initialise parameters
    f = a
    dx = 0.
    dy = 0.
    dz = 0.
    temp = np.zeros((sizex, sizey))
    index = 1.
    #% other loop parameters
    amaxvalue = a.max()
    pibyfour = np.pi/4.
    threetimespibyfour = 3.*pibyfour
    fivetimespibyfour = 5.*pibyfour
    seventimespibyfour = 7.*pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs((1./sinazimuth))
    dscos = np.abs((1./cosazimuth))
    tanaltitudebyscale = np.tan(altitude) / scale
    #% main loop
    while (amaxvalue >= dz and np.abs(dx) < sizex and np.abs(dy) < sizey):
        # if forsvf == 0:
        #     # dlg.progressBar.setValue(index)
    #while np.logical_and(np.logical_and(amaxvalue >= dz, np.abs(dx) <= sizex), np.abs(dy) <= sizey):(np.logical_and(amaxvalue >= dz, np.abs(dx) <= sizex), np.abs(dy) <= sizey):
        #if np.logical_or(np.logical_and(pibyfour <= azimuth, azimuth < threetimespibyfour), np.logical_and(fivetimespibyfour <= azimuth, azimuth < seventimespibyfour)):
        if (pibyfour <= azimuth and azimuth < threetimespibyfour or fivetimespibyfour <= azimuth and azimuth < seventimespibyfour):
            dy = signsinazimuth * index
            dx = -1. * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1. * signcosazimuth * index
            ds = dscos

        #% note: dx and dy represent absolute values while ds is an incremental value
        dz = ds *index * tanaltitudebyscale
        temp[0:sizex, 0:sizey] = 0.
        absdx = np.abs(dx)
        absdy = np.abs(dy)
        xc1 = (dx+absdx)/2.+1.
        xc2 = sizex+(dx-absdx)/2.
        yc1 = (dy+absdy)/2.+1.
        yc2 = sizey+(dy-absdy)/2.
        xp1 = -((dx-absdx)/2.)+1.
        xp2 = sizex-(dx+absdx)/2.
        yp1 = -((dy-absdy)/2.)+1.
        yp2 = sizey-(dy+absdy)/2.
        temp[int(xp1)-1:int(xp2), int(yp1)-1:int(yp2)] = a[int(xc1)-1:int(xc2), int(yc1)-1:int(yc2)]-dz
        # f = np.maximum(f, temp)  # bad performance in python3. Replaced with fmax
        f = np.fmax(f, temp)
        index += 1.

    f = f-a
    f = np.logical_not(f)
    sh = np.double(f)

    return sh

# @jit(nopython=True)


def shadowingfunction_20_updated(a, vegdem, vegdem2, azimuth, altitude, scale, amaxvalue, bush, forsvf):
    ''' Adapted by LB for standalone execution
    optimized for some redundant operations'''
    
    # Conversion
    degrees = np.pi / 180.
    azimuth = azimuth * degrees
    altitude = altitude * degrees
    
    # Measure the size of grid
    sizex, sizey = a.shape
    
    # Initialize parameters
    bushplant = bush > 1.
    sh = np.zeros((sizex, sizey))  # Shadows from buildings
    vbshvegsh = np.zeros((sizex, sizey))  # Vegetation blocking buildings
    vegsh = bushplant.astype(float)  # Vegetation shadow
    f = a.copy()
    
    pibyfour = np.pi / 4.
    threetimespibyfour = 3. * pibyfour
    fivetimespibyfour = 5. * pibyfour
    seventimespibyfour = 7. * pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs(1. / sinazimuth)
    dscos = np.abs(1. / cosazimuth)
    tanaltitudebyscale = np.tan(altitude) / scale
    index = 0
    
    # New case with pergola (thin vertical layer of vegetation), August 2021
    dzprev = 0

    # Main loop
    while amaxvalue >= dzprev and np.abs(dx) < sizex and np.abs(dy) < sizey:
        if (pibyfour <= azimuth < threetimespibyfour) or (fivetimespibyfour <= azimuth < seventimespibyfour):
            dy = signsinazimuth * index
            dx = -signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -signcosazimuth * index
            ds = dscos
        
        # Note: dx and dy represent absolute values while ds is an incremental value
        dz = (ds * index) * tanaltitudebyscale
        absdx, absdy = np.abs(dx), np.abs(dy)
        
        xc1, xc2 = int((dx + absdx) / 2.), int(sizex + (dx - absdx) / 2.)
        yc1, yc2 = int((dy + absdy) / 2.), int(sizey + (dy - absdy) / 2.)
        xp1, xp2 = int(-((dx - absdx) / 2.)), int(sizex - (dx + absdx) / 2.)
        yp1, yp2 = int(-((dy - absdy) / 2.)), int(sizey - (dy + absdy) / 2.)

        tempvegdem = vegdem[xc1:xc2, yc1:yc2] - dz
        tempvegdem2 = vegdem2[xc1:xc2, yc1:yc2] - dz
        temp = a[xc1:xc2, yc1:yc2] - dz

        f = np.fmax(f, temp)  # Moving building shadow
        sh = (f > a).astype(float)
        
        fabovea = tempvegdem > a[xp1:xp2, yp1:yp2]
        gabovea = tempvegdem2 > a[xp1:xp2, yp1:yp2]

        # New pergola condition
        templastfabovea = vegdem[xc1:xc2, yc1:yc2] - dzprev
        templastgabovea = vegdem2[xc1:xc2, yc1:yc2] - dzprev
        lastfabovea = templastfabovea > a[xp1:xp2, yp1:yp2]
        lastgabovea = templastgabovea > a[xp1:xp2, yp1:yp2]
        dzprev = dz

        vegsh2 = (fabovea | gabovea | lastfabovea | lastgabovea).astype(float)
        vegsh2[vegsh2 == 4] = 0.
        vegsh2[vegsh2 > 0] = 1.

        vegsh = np.fmax(vegsh, vegsh2)
        vegsh[(vegsh * sh > 0.)] = 0.
        vbshvegsh += vegsh

        index += 1

    sh = 1. - sh
    vbshvegsh[vbshvegsh > 0.] = 1.
    vbshvegsh -= vegsh
    vegsh = 1. - vegsh
    vbshvegsh = 1. - vbshvegsh

    shadowresult = {'sh': sh, 'vegsh': vegsh, 'vbshvegsh': vbshvegsh}

    return shadowresult

def shadowingfunction_20(a, vegdem, vegdem2, azimuth, altitude, scale, amaxvalue, bush, forsvf):
    ''' Adapted by LB for standalone execution'''
    # plt.ion()
    # fig = plt.figure(figsize=(24, 7))
    # plt.axis('image')
    # ax1 = plt.subplot(2, 3, 1)
    # ax2 = plt.subplot(2, 3, 2)
    # ax3 = plt.subplot(2, 3, 3)
    # ax4 = plt.subplot(2, 3, 4)
    # ax5 = plt.subplot(2, 3, 5)
    # ax6 = plt.subplot(2, 3, 6)
    # ax1.title.set_text('fabovea')
    # ax2.title.set_text('gabovea')
    # ax3.title.set_text('vegsh at ' + str(altitude))
    # ax4.title.set_text('lastfabovea')
    # ax5.title.set_text('lastgabovea')
    # ax6.title.set_text('vegdem')

    # This function casts shadows on buildings and vegetation units.
    # New capability to deal with pergolas 20210827

    # conversion
    degrees = np.pi/180.
    azimuth = azimuth * degrees
    altitude = altitude * degrees
    
    # measure the size of grid
    sizex = a.shape[0]
    sizey = a.shape[1]
    
    # progressbar for svf plugin
    if forsvf == 0:
        barstep = np.max([sizex, sizey])
        # dlg.progressBar.setRange(0, barstep)
        # dlg.progressBar.setValue(0)

    # initialise parameters
    dx = 0.
    dy = 0.
    dz = 0.
    temp = np.zeros((sizex, sizey))
    tempvegdem = np.zeros((sizex, sizey))
    tempvegdem2 = np.zeros((sizex, sizey))
    templastfabovea = np.zeros((sizex, sizey))
    templastgabovea = np.zeros((sizex, sizey))
    bushplant = bush > 1.
    sh = np.zeros((sizex, sizey)) #shadows from buildings
    vbshvegsh = np.zeros((sizex, sizey)) #vegetation blocking buildings
    vegsh = np.add(np.zeros((sizex, sizey)), bushplant, dtype=float) #vegetation shadow
    vegsh = np.zeros((sizex, sizey)) #vegetation shadow
    f = a

    pibyfour = np.pi / 4.
    threetimespibyfour = 3. * pibyfour
    fivetimespibyfour = 5.* pibyfour
    seventimespibyfour = 7. * pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs((1./sinazimuth))
    dscos = np.abs((1./cosazimuth))
    tanaltitudebyscale = np.tan(altitude) / scale
    # index = 1
    index = 0

    # new case with pergola (thin vertical layer of vegetation), August 2021
    dzprev = 0

    # main loop
    while (amaxvalue >= dz) and (np.abs(dx) < sizex) and (np.abs(dy) < sizey):
        # if forsvf == 0:
            # dlg.progressBar.setValue(index)
        if ((pibyfour <= azimuth) and (azimuth < threetimespibyfour) or (fivetimespibyfour <= azimuth) and (azimuth < seventimespibyfour)):
            dy = signsinazimuth * index
            dx = -1. * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1. * signcosazimuth * index
            ds = dscos
        # note: dx and dy represent absolute values while ds is an incremental value
        dz = (ds * index) * tanaltitudebyscale
        tempvegdem[0:sizex, 0:sizey] = 0.
        tempvegdem2[0:sizex, 0:sizey] = 0.
        temp[0:sizex, 0:sizey] = 0.
        templastfabovea[0:sizex, 0:sizey] = 0.
        templastgabovea[0:sizex, 0:sizey] = 0.
        absdx = np.abs(dx)
        absdy = np.abs(dy)
        xc1 = int((dx+absdx)/2.)
        xc2 = int(sizex+(dx-absdx)/2.)
        yc1 = int((dy+absdy)/2.)
        yc2 = int(sizey+(dy-absdy)/2.)
        xp1 = int(-((dx-absdx)/2.))
        xp2 = int(sizex-(dx+absdx)/2.)
        yp1 = int(-((dy-absdy)/2.))
        yp2 = int(sizey-(dy+absdy)/2.)


        tempvegdem[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dz
        tempvegdem2[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dz
        temp[xp1:xp2, yp1:yp2] = a[xc1:xc2, yc1:yc2]-dz

        f = np.fmax(f, temp) #Moving building shadow
        sh[(f > a)] = 1.
        sh[(f <= a)] = 0.
        fabovea = tempvegdem > a #vegdem above DEM
        gabovea = tempvegdem2 > a #vegdem2 above DEM
        
        #new pergola condition
        # templastfabovea[xp1:xp2, yp1:yp2] = tempvegdem[xc1:xc2, yc1:yc2]-dzprev
        # templastgabovea[xp1:xp2, yp1:yp2] = tempvegdem2[xc1:xc2, yc1:yc2]-dzprev
        templastfabovea[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2]-dzprev
        templastgabovea[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2]-dzprev
        lastfabovea = templastfabovea > a
        lastgabovea = templastgabovea > a
        dzprev = dz
        vegsh2 = np.add(np.add(np.add(fabovea, gabovea, dtype=float),lastfabovea, dtype=float),lastgabovea, dtype=float)
        vegsh2[vegsh2 == 4] = 0.
        # # vegsh2[vegsh2 == 1] = 0. # This one is the ultimate question...
        vegsh2[vegsh2 > 0] = 1.

        vegsh = np.fmax(vegsh, vegsh2)
        vegsh[(vegsh * sh > 0.)] = 0.
        vbshvegsh = vegsh + vbshvegsh # removing shadows 'behind' buildings

        # im1 = ax1.imshow(fabovea)
        # im2 = ax2.imshow(gabovea)
        # im3 = ax3.imshow(vegsh)
        # im4 = ax4.imshow(lastfabovea)
        # im5 = ax5.imshow(lastgabovea)
        # im6 = ax6.imshow(vegshtest)
        # im1 = ax1.imshow(tempvegdem)
        # im2 = ax2.imshow(tempvegdem2)
        # im3 = ax3.imshow(vegsh)
        # im4 = ax4.imshow(templastfabovea)
        # im5 = ax5.imshow(templastgabovea)
        # im6 = ax6.imshow(vegshtest)
        # plt.show()
        # plt.pause(0.05)

        index += 1.

    sh = 1.-sh
    vbshvegsh[(vbshvegsh > 0.)] = 1.
    vbshvegsh = vbshvegsh-vegsh
    vegsh = 1.-vegsh
    vbshvegsh = 1.-vbshvegsh

    # plt.close()
    # plt.ion()
    # fig = plt.figure(figsize=(24, 7))
    # plt.axis('image')
    # ax1 = plt.subplot(1, 3, 1)
    # im1 = ax1.imshow(vegsh)
    # plt.colorbar(im1)

    # ax2 = plt.subplot(1, 3, 2)
    # im2 = ax2.imshow(vegdem2)
    # plt.colorbar(im2)
    # plt.title('TDSM')

    # ax3 = plt.subplot(1, 3, 3)
    # im3 = ax3.imshow(vegdem)
    # plt.colorbar(im3)
    # plt.tight_layout()
    # plt.title('CDSM')
    # plt.show()
    # plt.pause(0.05)

    shadowresult = {'sh': sh, 'vegsh': vegsh, 'vbshvegsh': vbshvegsh}

    return shadowresult


def shadowingfunction_20_old(a, vegdem, vegdem2, azimuth, altitude, scale, amaxvalue, bush, dlg, forsvf):

    #% This function casts shadows on buildings and vegetation units
    #% conversion
    degrees = np.pi/180.
    if azimuth == 0.0:
        azimuth = 0.000000000001
    azimuth = np.dot(azimuth, degrees)
    altitude = np.dot(altitude, degrees)
    #% measure the size of the image
    sizex = a.shape[0]
    sizey = a.shape[1]
    #% initialise parameters
    if forsvf == 0:
        barstep = np.max([sizex, sizey])
        dlg.progressBar.setRange(0, barstep)
        dlg.progressBar.setValue(0)

    dx = 0.
    dy = 0.
    dz = 0.
    temp = np.zeros((sizex, sizey))
    tempvegdem = np.zeros((sizex, sizey))
    tempvegdem2 = np.zeros((sizex, sizey))
    sh = np.zeros((sizex, sizey))
    vbshvegsh = np.zeros((sizex, sizey))
    vegsh = np.zeros((sizex, sizey))
    tempbush = np.zeros((sizex, sizey))
    f = a
    g = np.zeros((sizex, sizey))
    bushplant = bush > 1.
    pibyfour = np.pi/4.
    threetimespibyfour = 3.*pibyfour
    fivetimespibyfour = 5.*pibyfour
    seventimespibyfour = 7.*pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs((1./sinazimuth))
    dscos = np.abs((1./cosazimuth))
    tanaltitudebyscale = np.tan(altitude) / scale
    index = 1

    #% main loop
    while (amaxvalue >= dz and np.abs(dx) < sizex and np.abs(dy) < sizey):
        if forsvf == 0:
            dlg.progressBar.setValue(index)
        if (pibyfour <= azimuth and azimuth < threetimespibyfour or fivetimespibyfour <= azimuth and azimuth < seventimespibyfour):
            dy = signsinazimuth * index
            dx = -1. * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1. * signcosazimuth * index
            ds = dscos
        #% note: dx and dy represent absolute values while ds is an incremental value
        dz = np.dot(np.dot(ds, index), tanaltitudebyscale)
        tempvegdem[0:sizex, 0:sizey] = 0.
        tempvegdem2[0:sizex, 0:sizey] = 0.
        temp[0:sizex, 0:sizey] = 0.
        absdx = np.abs(dx)
        absdy = np.abs(dy)
        xc1 = (dx+absdx)/2.+1.
        xc2 = sizex+(dx-absdx)/2.
        yc1 = (dy+absdy)/2.+1.
        yc2 = sizey+(dy-absdy)/2.
        xp1 = -((dx-absdx)/2.)+1.
        xp2 = sizex-(dx+absdx)/2.
        yp1 = -((dy-absdy)/2.)+1.
        yp2 = sizey-(dy+absdy)/2.
        tempvegdem[int(xp1)-1:int(xp2), int(yp1)-1:int(yp2)] = vegdem[int(xc1)-1:int(xc2), int(yc1)-1:int(yc2)]-dz
        tempvegdem2[int(xp1)-1:int(xp2), int(yp1)-1:int(yp2)] = vegdem2[int(xc1)-1:int(xc2), int(yc1)-1:int(yc2)]-dz
        temp[int(xp1)-1:int(xp2), int(yp1)-1:int(yp2)] = a[int(xc1)-1:int(xc2), int(yc1)-1:int(yc2)]-dz
        # f = np.maximum(f, temp) # bad performance in python3. Replaced with fmax
        f = np.fmax(f, temp)
        sh[(f > a)] = 1.
        sh[(f <= a)] = 0.
        #%Moving building shadow
        fabovea = tempvegdem > a
        #%vegdem above DEM
        gabovea = tempvegdem2 > a
        #%vegdem2 above DEM
        # vegsh2 = np.float(fabovea)-np.float(gabovea)
        vegsh2 = np.subtract(fabovea, gabovea, dtype=float)
        # vegsh = np.maximum(vegsh, vegsh2) # bad performance in python3. Replaced with fmax
        vegsh = np.fmax(vegsh, vegsh2)
        vegsh[(vegsh*sh > 0.)] = 0.
        #% removing shadows 'behind' buildings
        vbshvegsh = vegsh+vbshvegsh
        #% vegsh at high sun altitudes
        if index == 1.:
            firstvegdem = tempvegdem-temp
            firstvegdem[(firstvegdem <= 0.)] = 1000.
            vegsh[(firstvegdem < dz)] = 1.
            vegsh = vegsh*(vegdem2 > a)
            vbshvegsh = np.zeros((sizex, sizey))

        #% Bush shadow on bush plant
        if np.logical_and(bush.max() > 0., np.max((fabovea*bush)) > 0.):
            tempbush[0:sizex, 0:sizey] = 0.
            tempbush[int(xp1)-1:int(xp2), int(yp1)-1:int(yp2)] = bush[int(xc1)-1:int(xc2),int(yc1)-1:int(yc2)]-dz
            # g = np.maximum(g, tempbush) # bad performance in python3. Replaced with fmax
            g = np.fmax(g, tempbush)
            g *= bushplant
        index += 1.

    sh = 1.-sh
    vbshvegsh[(vbshvegsh > 0.)] = 1.
    vbshvegsh = vbshvegsh-vegsh

    if bush.max() > 0.:
        g = g-bush
        g[(g > 0.)] = 1.
        g[(g < 0.)] = 0.
        vegsh = vegsh-bushplant+g
        vegsh[(vegsh<0.)] = 0.

    vegsh[(vegsh > 0.)] = 1.
    vegsh = 1.-vegsh
    vbshvegsh = 1.-vbshvegsh

    shadowresult = {'sh': sh, 'vegsh': vegsh, 'vbshvegsh': vbshvegsh}

    return shadowresult