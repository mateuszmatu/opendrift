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

o.set_config('general:deac', True)



o.set_config('drift:vertical_mixing', False)
o.set_config('drift:vertical_advection', False)
max = 999
min = -999

o.set_config('deac:max', max)
o.set_config('deac:min', min)
o.set_config('deac:method', 'hard_minmax')
o.set_config('deac:variable', 'light')

time = datetime(2025, 10, 25, 12)

pos = [7.3, 57.2]
o.seed_elements(pos[0], pos[1], z=-15, number=1,
                time=time, prefered_light=10, vertical_swim_speed=0.001)
o.run(duration=timedelta(hours=10), time_step=timedelta(minutes=15), time_step_output=timedelta(minutes=15), outfile='../light_attraction.nc')                                   
