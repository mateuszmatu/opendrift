#!/usr/bin/env python
"""
Cod egg
=============
"""

from opendrift.readers import reader_netCDF_CF_generic
from opendrift.models.pelagicegghealth import PelagicEggHealthDrift
from datetime import datetime, timedelta

o = PelagicEggHealthDrift(loglevel=20)  # Set loglevel to 0 for debug information

# Forcing with Topaz ocean model and MEPS atmospheric model
#o.add_readers_from_list([
#    'https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be',
#    'https://thredds.met.no/thredds/dodsC/mepslatest/meps_lagged_6_h_latest_2_5km_latest.nc'])
#

o.add_reader(reader_netCDF_CF_generic.Reader('https://thredds.met.no/thredds/dodsC/fou-hi/norkystv3_800m_m00_be'))
#%%
# Adjusting some configuration
o.set_config('drift:vertical_mixing', False)
#o.set_config('vertical_mixing:diffusivitymodel', 'windspeed_Sundby1983') # windspeed parameterization for eddy diffusivity

#o.set_config('general:premature_deactivation', 'exposure')
#o.set_config('general:deactivation_exposure', 'simple_minmax_exposure')

o.set_config('general:premature_deactivation', 'shortwave_minmax')


max = 999
min = -999
o.set_config('general:deactivation_max', max)
o.set_config('general:deactivation_min', min)

#%%
# Vertical mixing requires fast time step
o.set_config('vertical_mixing:timestep', 60.) # seconds

#%%
# spawn NEA cod eggs at defined position and time
time = datetime(2025, 9, 25, 12)

pos = [7.3, 57.2]
o.seed_elements(pos[0], pos[1], z=-15, radius=100, number=2,
                time=time, diameter=0.0014, neutral_buoyancy_salinity=31.25)

#%%
# Running model
o.run(duration=timedelta(hours=5), time_step=3600, outfile='../test.nc')                                   

#%%
# Print and plot results.
# At the end the wind vanishes, and eggs come to surface
#print(o)

#o.plot(fast=True)
#o.animation(fast=True, color='z')

#%%
# .. image:: /gallery/animations/example_codegg_0.gif

#%% Vertical distribution of particles
#o.animate_vertical_distribution()

#%%
# .. image:: /gallery/animations/example_codegg_1.gif
