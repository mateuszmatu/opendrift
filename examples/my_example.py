from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.my_element import MyElementDrift
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import xarray as xr

o = MyElementDrift(loglevel=20)  


o.add_reader(reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'))

#o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity

#o.set_config('general:premature_deactivation', 'exposure')
#o.set_config('general:deactivation_exposure', 'simple_minmax_exposure')

o.set_config('general:deac', 'shortwave_minmax')



o.set_config('drift:vertical_mixing', True)
max = 999
min = -999

o.set_config('general:deac_max', max)
o.set_config('general:deac_min', min)
o.set_config('vertical_mixing:timestep', 60.) # seconds

time = datetime(2025, 10, 25, 12)

pos = [7.3, 57.2]
o.seed_elements(pos[0], pos[1], z=-15, radius=100, number=2,
                time=time)
o.run(duration=timedelta(hours=5), time_step=3600, outfile='../test.nc')                                   

ds = xr.open_dataset('../test.nc')
plt.scatter(ds.lon, ds.lat)
plt.show()