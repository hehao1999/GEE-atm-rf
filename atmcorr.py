import ee
from Py6S import *
import datetime
import math
import os
import sys
from cloudmask import maskS2clouds
sys.path.append(os.path.join(os.path.dirname(os.getcwd()),'bin'))
from atmospheric import Atmospheric


"""初始化GEE"""
os.environ['HTTP_PROXY'] = 'http://127.0.0.1:10080'
os.environ['HTTPS_PROXY'] = 'https://127.0.0.1:10080'
# ee.Authenticate() # 第一次登录时使用
ee.Initialize()
print(1)
"""导入研究区"""
area = ee.FeatureCollection("users/hehao/area").geometry()
imgs = ee.ImageCollection('COPERNICUS/S2')\
    .filterDate('2019-08-01', '2019-08-31')\
    .filterBounds(area)\
    .sort('CLOUDY_PIXEL_PERCENTAGE', False)\
    .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 40))

# **The first Sentinel 2 image
S2 = ee.Image("COPERNICUS/S2/20190827T023551_20190827T024252_T51TXH")


# metadata
info = S2.getInfo()['properties']
scene_date = datetime.datetime.utcfromtimestamp(info['system:time_start']/1000)# i.e. Python uses seconds, EE uses milliseconds
date = ee.Date(scene_date)
solar_z = info['MEAN_SOLAR_ZENITH_ANGLE']

#需要修改
geom = ee.Geometry.Point(125.266, 42.86)

toa = S2.divide(10000)

print(2)

"""atmospheric constituents"""
h2o = Atmospheric.water(geom,date).getInfo()
o3 = Atmospheric.ozone(geom,date).getInfo()
aot = Atmospheric.aerosol(geom,date).getInfo()

"""target altitude (km)"""
SRTM = ee.Image('CGIAR/SRTM90_V4')# Shuttle Radar Topography mission covers *most* of the Earth
alt = SRTM.reduceRegion(reducer = ee.Reducer.mean(),geometry = geom.centroid()).get('elevation').getInfo()
km = alt/1000 # i.e. Py6S uses units of kilometers

"""6S object"""
# Instantiate
s = SixS()

# Atmospheric constituents
s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o,o3)
s.aero_profile = AeroProfile.Continental
s.aot550 = aot

# Earth-Sun-satellite geometry
s.geometry = Geometry.User()
s.geometry.solar_z = solar_z        # solar zenith angle
s.geometry.month = scene_date.month # month and day used for Earth-Sun distance
s.geometry.day = scene_date.day     # month and day used for Earth-Sun distance
s.altitudes.set_sensor_satellite_level()
s.altitudes.set_target_custom_altitude(km)

"""Spectral Response functions"""


def spectralResponseFunction(bandname):
    """
    Extract spectral response function for given band name
    """

    bandSelect = {
        'B1': PredefinedWavelengths.S2A_MSI_01,
        'B2': PredefinedWavelengths.S2A_MSI_02,
        'B3': PredefinedWavelengths.S2A_MSI_03,
        'B4': PredefinedWavelengths.S2A_MSI_04,
        'B5': PredefinedWavelengths.S2A_MSI_05,
        'B6': PredefinedWavelengths.S2A_MSI_06,
        'B7': PredefinedWavelengths.S2A_MSI_07,
        'B8': PredefinedWavelengths.S2A_MSI_08,
        'B8A': PredefinedWavelengths.S2A_MSI_09,
        'B9': PredefinedWavelengths.S2A_MSI_10,
        'B10': PredefinedWavelengths.S2A_MSI_11,
        'B11': PredefinedWavelengths.S2A_MSI_12,
        'B12': PredefinedWavelengths.S2A_MSI_13,
    }

    return Wavelength(bandSelect[bandname])

"""TOA Reflectance to Radiance"""


def toa_to_rad(bandname):
    """
    Converts top of atmosphere reflectance to at-sensor radiance
    """

    # solar exoatmospheric spectral irradiance
    ESUN = info['SOLAR_IRRADIANCE_' + bandname]
    solar_angle_correction = math.cos(math.radians(solar_z))

    # Earth-Sun distance (from day of year)
    doy = scene_date.timetuple().tm_yday
    d = 1 - 0.01672 * math.cos(0.9856 * (
                doy - 4))  # http://physics.stackexchange.com/questions/177949/earth-sun-distance-on-a-given-day-of-the-year

    # conversion factor
    multiplier = ESUN * solar_angle_correction / (math.pi * d ** 2)

    # at-sensor radiance
    rad = toa.select(bandname).multiply(multiplier)

    return rad

"""Radiance to Surface Reflectance"""


def surface_reflectance(bandname):
    """
    Calculate surface reflectance from at-sensor radiance given waveband name
    """

    # **Set remaining geometry parameters
    s.geometry.view_a = info['MEAN_INCIDENCE_AZIMUTH_ANGLE_' + bandname]
    s.geometry.view_z = info['MEAN_INCIDENCE_ZENITH_ANGLE_' + bandname]

    # run 6S for this waveband
    s.wavelength = spectralResponseFunction(bandname)
    s.run()

    # extract 6S outputs
    Edir = s.outputs.direct_solar_irradiance  # direct solar irradiance
    Edif = s.outputs.diffuse_solar_irradiance  # diffuse solar irradiance
    Lp = s.outputs.atmospheric_intrinsic_radiance  # path radiance
    absorb = s.outputs.trans['global_gas'].upward  # absorption transmissivity
    scatter = s.outputs.trans['total_scattering'].upward  # scattering transmissivity
    tau2 = absorb * scatter  # total transmissivity

    # radiance to surface reflectance
    rad = toa_to_rad(bandname)
    ref = rad.subtract(Lp).multiply(math.pi).divide(tau2 * (Edir + Edif))

    return ref

"""Atmospheric Correction"""
# surface reflectance rgb
b = surface_reflectance('B2')
g = surface_reflectance('B3')
r = surface_reflectance('B4')
ref = r.addBands(g).addBands(b)

# all wavebands
output = S2.select('QA60')
for band in ['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12']:
    output = output.addBands(surface_reflectance(band))

"""Display results"""
from IPython.display import display, Image

region = geom.buffer(5000).bounds().getInfo()['coordinates']
channels = ['B4','B3','B2']

original = Image(url=toa.select(channels).getThumbUrl({
                'region':region,
                'min':0,
                'max':0.25
                }))

corrected = Image(url=ref.select(channels).getThumbUrl({
                'region':region,
                'min':0,
                'max':0.25
                }))

display(original, corrected)

"""Export to Asset"""
# set some properties for export
dateString = scene_date.strftime("%Y-%m-%d")
output = output.set({'satellite':'Sentinel 2',
              'fileID':info['system:index'],
              'date':dateString,
              'aerosol_optical_thickness':aot,
              'water_vapour':h2o,
              'ozone':o3})

# **define my assetID
# **increment the number for each image processed
assetID = 'users/hehao/' + 'Sentinel_Corrected_' + 'ttt'

# export
export = ee.batch.Export.image.toAsset(
    image=output,
    description='sentinel2_atmcorr_export',
    #region = region,
    assetId = assetID,
    maxPixels = 10000000000000
)

export.start()
print('OK')