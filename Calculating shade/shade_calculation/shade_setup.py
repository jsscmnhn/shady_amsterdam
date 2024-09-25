import os
import datetime as dt
import sun_position as sp
import shadowingfunctions as shadow
from osgeo import gdal, osr
from osgeo.gdalconst import *
import numpy as np
import geopandas as gpd
import pandas as pd



### Shade calculation setup

def shadecalculation_setup(filepath_dsm='None', filepath_veg='None',
                           date=dt.datetime.now(),
                           intervalTime=30, onetime=1, filepath_save='None', UTC=0, dst=1, useveg=0, trunkheight=25,
                           transmissivity=20):
    '''Calculates spot, hourly and or daily shading for a DSM
    Needs:
    filepath_dsm = a path to a (building) dsm,
    date = a datetime object (defaults to datetime.now()),
    intervaltime = an integer in minutes as interval,
    onetime = to differentiate between single shade cast and day (default is single timestamp)
    filepath_save = 'path to the folder in which to save the files'
    UTC = defaults to 0, optional in plugin, denotes the UCT shift (AMS +2)
    dst = defaults to 1, not sure what this does. 1 or 0
    '''
    print(date)

    gdal_dsm = gdal.Open(filepath_dsm)
    dsm = gdal_dsm.ReadAsArray().astype(float)

    # response to issue #85
    nd = gdal_dsm.GetRasterBand(1).GetNoDataValue()
    dsm[dsm == nd] = 0.
    if dsm.min() < 0:
        dsm = dsm + np.abs(dsm.min())

    sizex = dsm.shape[0]
    sizey = dsm.shape[1]

    old_cs = osr.SpatialReference()
    dsm_ref = gdal_dsm.GetProjection()
    print(dsm_ref)
    # dsm_ref = dsmlayer.crs().toWkt()
    old_cs.ImportFromWkt(dsm_ref)

    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""

    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(wgs84_wkt)

    transform = osr.CoordinateTransformation(old_cs, new_cs)

    width = gdal_dsm.RasterXSize
    height = gdal_dsm.RasterYSize
    gt = gdal_dsm.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width * gt[4] + height * gt[5]
    lonlat = transform.TransformPoint(minx, miny)
    geotransform = gdal_dsm.GetGeoTransform()
    scale = 1 / geotransform[1]

    gdalver = float(gdal.__version__[0])
    if gdalver >= 3.:
        lon = lonlat[1]  # changed to gdal 3
        lat = lonlat[0]  # changed to gdal 3
    else:
        lon = lonlat[0]  # changed to gdal 2
        lat = lonlat[1]  # changed to gdal 2
    # print('lon:' + str(lon))
    # print('lat:' + str(lat))

    ## Import vegetation dsm
    trans = transmissivity / 100.0
    if useveg == 1:
        usevegdem = 1

        # Changed import conditions
        if filepath_veg == "None":
            print("No vegetation filepath given")
            return

        dataSet = gdal.Open(filepath_veg)
        vegdsm = dataSet.ReadAsArray().astype(float)

        vegsizex = vegdsm.shape[0]
        vegsizey = vegdsm.shape[1]

        if not (vegsizex == sizex) & (vegsizey == sizey):
            print("Error; All grids must be of same extent and resolution")
            return
        ## Defaulted trunk zone to 25%
        # else:
        trunkratio = trunkheight / 100.0
        vegdsm2 = vegdsm * trunkratio

        vegsizex = vegdsm2.shape[0]
        vegsizey = vegdsm2.shape[1]

        if not (vegsizex == sizex) & (vegsizey == sizey):  # &
            print("Error; All grids must be of same extent and resolution")
            return
    else:
        vegdsm = 0
        vegdsm2 = 0
        usevegdem = 0

    wallsh = 0
    wheight = 0
    waspect = 0

    if filepath_save == 'None':
        print("Error", "No output path given")
        return
    else:
        # date = self.dlg.calendarWidget.selectedDate()
        year = date.year
        month = date.month
        day = date.day
        # TODO: Check what UTC does. Optional in plugin, so set to 0
        # UTC = self.dlg.spinBoxUTC.value()
        if onetime == 1:
            time = date.time()
            print(time)
            hour = time.hour
            minu = time.minute
            sec = time.second
        else:
            onetime = 0
            hour = 0
            minu = 0
            sec = 0

        tv = [year, month, day, hour, minu, sec]
        # TODO: What exactly is set by timeInterval. For now this just takes the minutes as set above.
        # intervalTime = self.dlg.intervalTimeEdit.time()
        # self.timeInterval = intervalTime.minute() + (intervalTime.hour() * 60) + (intervalTime.second()/60)

        shadowresult = dailyshading(dsm, vegdsm, vegdsm2, scale, lon, lat, sizex, sizey, tv, UTC, usevegdem,
                                    intervalTime, onetime, filepath_save, gdal_dsm, trans,
                                    dst, wallsh, wheight, waspect)

        shfinal = shadowresult["shfinal"]
        time_vector = shadowresult["time_vector"]

        if onetime == 0:
            timestr = time_vector.strftime("%Y%m%d")
            # Changed filepath to include the tile_id in the filename.
            # savestr = '/shadow_fraction_on_'
            # savestr = '/shadow_fraction_on_'
            savestr = 'shadow_'
        else:
            timestr = time_vector.strftime("%Y%m%d_%H%M")
            # savestr = '/Shadow_at_'
            savestr = 'shadow_'

    filename = filepath_save + savestr + timestr + '.tif'

    ## TODO: change to saverasternd or other function
    saveraster(gdal_dsm, filename, shfinal)


############## DAILYSHADING ################

def dailyshading(dsm, vegdsm, vegdsm2, scale, lon, lat, sizex, sizey, tv, UTC, usevegdem, timeInterval, onetime, folder,
                 gdal_data, trans, dst, wallshadow, wheight, waspect):
    # lon = lonlat[0]
    # lat = lonlat[1]
    year = tv[0]
    month = tv[1]
    day = tv[2]

    alt = np.median(dsm)
    location = {'longitude': lon, 'latitude': lat, 'altitude': alt}

    if usevegdem == 1:
        psi = trans
        # amaxvalue
        vegmax = vegdsm.max()
        amaxvalue = dsm.max() - dsm.min()
        amaxvalue = np.maximum(amaxvalue, vegmax)

        # Elevation vegdsms if buildingDSM includes ground heights
        vegdem = vegdsm + dsm
        vegdem[vegdem == dsm] = 0
        vegdem2 = vegdsm2 + dsm
        vegdem2[vegdem2 == dsm] = 0

        # Bush separation
        bush = np.logical_not((vegdem2 * vegdem)) * vegdem


    shtot = np.zeros((sizex, sizey))

    if onetime == 1:
        itera = 1
    else:
        itera = int(np.round(1440 / timeInterval))

    alt = np.zeros(itera)
    azi = np.zeros(itera)
    hour = int(0)
    index = 0
    time = dict()
    time['UTC'] = UTC


    walls = np.zeros((sizex, sizey))
    dirwalls = np.zeros((sizex, sizey))

    for i in range(0, itera):
        if onetime == 0:
            minu = int(timeInterval * i)
            if minu >= 60:
                hour = int(np.floor(minu / 60))
                minu = int(minu - hour * 60)
        else:
            minu = tv[4]
            hour = tv[3]

        doy = day_of_year(year, month, day)

        ut_time = doy - 1. + ((hour - dst) / 24.0) + (minu / (60. * 24.0)) + (0. / (60. * 60. * 24.0))

        if ut_time < 0:
            year = year - 1
            month = 12
            day = 31
            doy = day_of_year(year, month, day)
            ut_time = ut_time + doy - 1

        HHMMSS = dectime_to_timevec(ut_time)

        time['year'] = year
        time['month'] = month
        time['day'] = day
        time['hour'] = HHMMSS[0]
        time['min'] = HHMMSS[1]
        time['sec'] = HHMMSS[2]

        sun = sp.sun_position(time, location)
        alt[i] = 90. - sun['zenith']
        azi[i] = sun['azimuth']

        if time['sec'] == 59:  # issue 228 and 256
            time['sec'] = 0
            time['min'] = time['min'] + 1
            if time['min'] == 60:
                time['min'] = 0
                time['hour'] = time['hour'] + 1
                if time['hour'] == 24:
                    time['hour'] = 0

        time_vector = dt.datetime(year, month, day, time['hour'], time['min'], time['sec'])
        timestr = time_vector.strftime("%Y%m%d_%H%M")

        if alt[i] > 0:
            if usevegdem == 0:
                sh = shadow.shadowingfunctionglobalradiation(dsm, azi[i], alt[i], scale, 0)
                # shtot = shtot + sh
            else:
                # changed to "optimized" function
                shadowresult = shadow.shadowingfunction_20(dsm, vegdem, vegdem2, azi[i], alt[i], scale, amaxvalue,
                                                           bush, 0)

                vegsh = shadowresult["vegsh"]
                sh = shadowresult["sh"]
                sh = sh - (1 - vegsh) * (1 - psi)

            if onetime == 0:
                filename = folder + '_Shadow_' + timestr + '_LST.tif'
                ## EDITED
                saveraster(gdal_data, filename, sh)

        shtot = shtot + sh
        index += 1

    shfinal = shtot / index
    shadowresult = {'shfinal': shfinal, 'time_vector': time_vector}

    # dlg.progressBar.setValue(0)

    return shadowresult


def day_of_year(yy, month, day):
    if (yy % 4) == 0:
        if (yy % 100) == 0:
            if (yy % 400) == 0:
                leapyear = 1
            else:
                leapyear = 0
        else:
            leapyear = 1
    else:
        leapyear = 0

    if leapyear == 1:
        dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    else:
        dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    doy = np.sum(dayspermonth[0:month - 1]) + day

    return doy


def dectime_to_timevec(dectime):
    # This subroutine converts dectime to individual hours, minutes and seconds

    doy = np.floor(dectime)

    DH = dectime - doy
    HOURS = int(24 * DH)

    DM = 24 * DH - HOURS
    MINS = int(60 * DM)

    DS = 60 * DM - MINS
    SECS = int(60 * DS)

    return (HOURS, MINS, SECS)


# def saveraster(gdal_data, filename, raster):
#     rows = gdal_data.RasterYSize
#     cols = gdal_data.RasterXSize
#
#     outDs = gdal.GetDriverByName("GTiff").Create(filename, cols, rows, int(1), GDT_Float32)
#     outBand = outDs.GetRasterBand(1)
#
#     # write the data
#     outBand.WriteArray(raster, 0, 0)
#     # flush data to disk, set the NoData value and calculate stats
#     outBand.FlushCache()
#     outBand.SetNoDataValue(-9999)
#
#     # georeference the image and set the projection
#     outDs.SetGeoTransform(gdal_data.GetGeoTransform())
#     outDs.SetProjection(gdal_data.GetProjection())

def saveraster(gdal_data, filename, raster):
    rows = gdal_data.RasterYSize
    cols = gdal_data.RasterXSize

    # Check if rows and cols are valid
    if rows <= 0 or cols <= 0:
        print(f"Invalid raster dimensions: rows={rows}, cols={cols}")
        return

    # Create the output dataset
    driver = gdal.GetDriverByName("GTiff")
    if driver is None:
        print("GTiff driver is not available.")
        return

    outDs = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)  # 1 = single band

    # Check if outDs was created successfully
    if outDs is None:
        print(f"Failed to create raster dataset for file: {filename}")
        return

    outBand = outDs.GetRasterBand(1)
    if outBand is None:
        print(f"Failed to get raster band for file: {filename}")
        return

    # Write the data
    try:
        outBand.WriteArray(raster)
    except Exception as e:
        print(f"Error writing raster data: {e}")
        return

    # Flush cache and set NoData value
    outBand.FlushCache()
    outBand.SetNoDataValue(-9999)

    # Georeference the image and set projection
    outDs.SetGeoTransform(gdal_data.GetGeoTransform())
    outDs.SetProjection(gdal_data.GetProjection())

    # Close the dataset properly to ensure the data is written to disk
    outDs = None
    print(f"Raster saved successfully to {filename}")

