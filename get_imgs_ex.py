import ee
from Py6S import *
import datetime
import math
import os
import sys
from atmospheric import Atmospheric

sys.path.append(os.path.join(os.path.dirname(os.getcwd()), 'bin'))

"""初始化GEE"""
os.environ['HTTP_PROXY'] = 'http://127.0.0.1:10080'
os.environ['HTTPS_PROXY'] = 'https://127.0.0.1:10080'
#ee.Authenticate() # 第一次登录时使用
ee.Initialize()
print("GEE Initialize Success")

"""导入研究区"""
area = ee.FeatureCollection("users/hehao/area").geometry() # 研究区矢量边界

# 筛选出研究区内符合条件的S2影像集
images = [
    ee.Image('COPERNICUS/S2/20190803T034539_20190803T035820_T48STD'),
    ee.Image('COPERNICUS/S2/20190811T035541_20190811T040503_T48STD'),
    ee.Image('COPERNICUS/S2/20190803T034539_20190803T035820_T48SUD'),
    ee.Image('COPERNICUS/S2/20190803T034539_20190803T034546_T48STE'),
    ee.Image('COPERNICUS/S2/20190811T035541_20190811T040503_T48STF'),
    ee.Image('COPERNICUS/S2/20190801T035541_20190801T035541_T48SUG'),
    ee.Image('COPERNICUS/S2/20190810T033539_20190810T034635_T48SWF'),
    ee.Image('COPERNICUS/S2/20190810T033539_20190810T034635_T48SXG'),
    ee.Image('COPERNICUS/S2/20190810T033539_20190810T034635_T49SBB'),
    ee.Image('COPERNICUS/S2/20190807T032549_20190807T033637_T49SCC'),
    ee.Image('COPERNICUS/S2/20190827T032539_20190827T033131_T49TEE'),
    ee.Image('COPERNICUS/S2/20190814T031549_20190814T031915_T49TFE'),
    ee.Image('COPERNICUS/S2/20190814T031549_20190814T031915_T49SFD'),
    #    ee.Image('COPERNICUS/S2/20190803T025551_20190803T030454_T50TPP'),
    ee.Image('COPERNICUS/S2/20190815T024549_20190815T025048_T50TQM'),
    ee.Image('COPERNICUS/S2/20190822T023559_20190822T023554_T51TVG')]
col = ee.ImageCollection.fromImages(images)

"""获取影像id列表"""
col_info = col.getInfo()['features']
img_id_list = []
for img_info in col_info:
    img_id_list.append(img_info['properties']['system:index'])

atm_img_list = []  # 用于存放进行过大气校正的影像

# 对集合内影像进行大气校正，由jupyter notebook内的代码改编为批处理
for i, img_id in enumerate(img_id_list):

    """metadata"""
    S2 = ee.Image('COPERNICUS/S2/' + img_id)
    info = S2.getInfo()['properties']
    scene_date = datetime.datetime.utcfromtimestamp(info['system:time_start'] / 1000)
    date = ee.Date(scene_date)
    solar_z = info['MEAN_SOLAR_ZENITH_ANGLE']
    coordinates = info['system:footprint']['coordinates']
    lat = []
    lon = []
    # 这里以影像最大外接矩形中心作为校正基点，获取该点的水气情况
    for coordinate in coordinates:
        lat.append(coordinate[0])
        lon.append(coordinate[1])
    geom = ee.Geometry.Point((max(lat) + min(lat)) / 2, (max(lon) + min(lon)) / 2)
    toa = S2.divide(10000)

    """atmospheric constituents"""
    h2o = Atmospheric.water(geom, date).getInfo()
    o3 = Atmospheric.ozone(geom, date).getInfo()
    aot = Atmospheric.aerosol(geom, date).getInfo()

    """target altitude (km)"""
    SRTM = ee.Image('CGIAR/SRTM90_V4')  # Shuttle Radar Topography mission covers *most* of the Earth
    alt = SRTM.reduceRegion(reducer=ee.Reducer.mean(), geometry=geom.centroid()).get('elevation').getInfo()
    km = alt / 1000  # i.e. Py6S uses units of kilometers

    """6S object"""
    # Instantiate
    s = SixS()
    # Atmospheric constituents
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o, o3)
    s.aero_profile = AeroProfile.Continental
    s.aot550 = aot
    # Earth-Sun-satellite geometry
    s.geometry = Geometry.User()
    s.geometry.solar_z = solar_z  # solar zenith angle
    s.geometry.month = scene_date.month  # month and day used for Earth-Sun distance
    s.geometry.day = scene_date.day  # month and day used for Earth-Sun distance
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
    # 保留还需要的波段
    output = surface_reflectance('B1')  # S2.select('QA60') # 也可以导入QA60波段，如果先进行去云处理需要修改后面的裁剪拼接代码
    for band in ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']:
        output = output.addBands(surface_reflectance(band))
    atm_img_list.append(output)
    print(img_id + '：success-' + str(i))
print('atmcorr success')


# 将图像转换为DOUBLE型，否则导出会发生flaot类型范围错误，不知为何只能在图像拼接以前进行toDouble处理
def to_double(image):
    return image.toDouble()


"""mosaic and clip"""
atm_col = ee.ImageCollection.fromImages(atm_img_list)
image = atm_col \
    .sort('CLOUDY_PIXEL_PERCENTAGE', False) \
    .map(to_double)\
    .mosaic() \
    .clip(area)


"""Export to Asset"""


def export_img(bandname, scale):
    export = ee.batch.Export.image.toAsset(
        image=image.select(bandname),
        description='extra' + bandname,     # task任务名称
        assetId='users/hehao/' + 'extra' + bandname,    # 导出影像名称
        scale=scale,
        maxPixels=1e13
    )
    export.start()
    print('export {0}'.format(bandname))

# 由于各个波段分辨率不同，分波段导出
# 未能将影像合成后导出为一幅多波段影像，因为分辨率不同，会导出为设置的分辨率
for band in [['B2',10], ['B3',10], ['B4',10], ['B5',20], ['B6',20],
             ['B7',20], ['B8',10], ['B8A',20], ['B11',20], ['B12',20]]:
    export_img(band[0], band[1])

print('All is OK')

