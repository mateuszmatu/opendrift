import numpy as np
import logging; logger = logging.getLogger(__name__)

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

class MyElement(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs
    """

    variables = Lagrangian3DArray.add_variables([
        ('health', {'dtype': np.float32,
                     'units': '',
                     'default': 100.}),
        ('light', {'dtype': np.float32,
                   'units': 'W m^-2',
                   'default': 0}),
        ('prefered_light', {'dtype': np.float32,
                           'units': 'umol m^-2 s^-1',
                           'default': 0}),
        ('vertical_swim_speed', {'dtype': np.float32,
                                  'units': 'm s^-1',
                                  'default': 0.0005}), # From Sandvik et. al. 2020
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # for NEA Co
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': '[]',
                                       'default': 37.25}),  # for NEA Cod
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 1028.})
    ])

class MyElementDrift(OceanDrift):

    ElementType = MyElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0,
            'store_previous_if': ['drift:vertical_advection', 'is', True]},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'upward_sea_water_velocity': {'fallback': 0,
            'skip_if': ['drift:vertical_advection', 'is', False]},
        'ocean_vertical_diffusivity': {'fallback': 0,
             'skip_if': ['drift:vertical_mixing', 'is', False],
             'profiles': True},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0,
            'skip_if': ['drift:stokes_drift', 'is', False]},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0,
            'skip_if': ['drift:stokes_drift', 'is', False]},
        'ocean_mixed_layer_thickness': {
            'fallback': 50, 'skip_if': ['drift:vertical_mixing', 'is', False]},
        'sea_floor_depth_below_sea_level': {'fallback': 10000},
        'land_binary_mask': {'fallback': None},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'net_downward_shortwave_flux_at_sea_water_surface': {'fallback': 150},
      }

    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta'}

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(MyElementDrift, self).__init__(*args, **kwargs)

        # By default, eggs do not strand towards coastline
        self._set_config_default('general:coastline_action', 'previous')

        # Vertical mixing is enabled by default, also at surface. Also vertical advection at surface.
        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:vertical_mixing_at_surface', True)
        self._set_config_default('drift:vertical_advection_at_surface', True)

        self._add_config({
            'general:deac': {'type': 'bool', 
                             'default': False,
                             'description': 'Turn on deactivation functions',
                             'level': CONFIG_LEVEL_ADVANCED},
            'deac:min':     {'type': 'float',
                            'default': -999.0, 'min': -999, 'max': 999, 
                            'units': 'None',
                            'description': 'Min threshold for deactivation function', 
                            'level': CONFIG_LEVEL_ADVANCED},
            'deac:max':     {'type': 'float',
                            'default': 999.0, 'min': -999, 'max': 999, 
                            'units': 'None',
                            'description': 'Max threshold for deactivation function', 
                            'level': CONFIG_LEVEL_ADVANCED},
            'deac:variable': {'type': 'str',
                              'default': 'sea_water_temperature',
                              'min_length': 1,
                              'max_length': 999,
                              'description': 'Physical variable for deactivation function', 
                              'level': CONFIG_LEVEL_ADVANCED},
            'deac:method':  {'type': 'enum',
                            'default': 'hard_minmax',
                            'description': 'Deactivation method',
                            'enum': ['hard_minmax', 'exposure'],
                            'level': CONFIG_LEVEL_ADVANCED},
            'deac:health_drain': {'type': 'float',
                                  'units': 'None',
                                  'default': 100,
                                  'description': 'Value to subtract from 100 each time step particle is outside of threshold. Deactivates particle at 0.',
                                  'min': 0,
                                  'max': 100,
                                  'level': CONFIG_LEVEL_ADVANCED},
            'my_element:avoided_salinity': {'type':'float', 'default':28,
                        'min':0, 'max':50, 'units': 'PSU',
                        'description': 'Salinity actively avoided',
                        'level': CONFIG_LEVEL_BASIC}
            })

        if self.get_config('deac:variable') not in self.required_variables and self.get_config('deac:variable') not in self.elements.variables.keys():
            raise ValueError(f'Variable {self.get_config('deac:variable')} is not in list of required variables or element variables.\n Add it with "OceanDrift.required_variables.update".')

    def deac(self):
        # need to rework this later
        deac_indices = []
        
        if self.get_config('deac:variable') in self.required_variables:
            considered_value = self.required_variables[self.get_config('deac:variable')]
        elif self.get_config('deac:variable') in self.elements.variables.keys():
            considered_value=self.elements.light #TODO get this to work to be more generic
            #considered_value = self.elements[self.get_config('deac:variable')]
        health_indices = [el < self.get_config('deac:min') 
                            or el > self.get_config('deac:max')
                            for el in considered_value]
        self.elements.health[health_indices == np.True_] -= self.get_config('deac:health_drain')

        deac_indices = [el <= 0 for el in self.elements.health]            

        if len(deac_indices) > 0:
            self.deactivate_elements(deac_indices, 'Deactivated.')

    def light_along_trajectory(self):
        self.elements.light = self.shortwave_radiation_at_depth(self.environment.net_downward_shortwave_flux_at_sea_water_surface, self.elements.z)

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        W = self.velocity_light()
        #W = self.velocity_shape()
        self.elements.terminal_velocity = W


    def detect_salinity(self):
        # Wants to be as close to light preference as possible without crossing salinity boundary.

        # Need salinity profile
        # Currently only senses the layer directly above or below. 

        Sprofiles = self.environment_profiles['sea_water_salinity']
        Tprofiles = self.environment_profiles['sea_water_temperature']

        

    def velocity_light(self):
        #µmol/m²/s = W/m² * 4.6

        photon = self.elements.light * 4.6
        W = np.where(self.elements.light == self.elements.prefered_light, 0, np.where(self.elements.light > self.elements.prefered_light, -self.elements.vertical_swim_speed, self.elements.vertical_swim_speed))
        
        return W

    def interpolate_profiles(self):
        Sprofiles = self.environment_profiles['sea_water_salinity']
        Tprofiles = self.environment_profiles['sea_water_temperature']

        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.uint8), 0)
            lower = np.minimum(upper+1, Tprofiles.shape[0]-1)
            weight_upper = 1 - (zi - upper)
        
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
                weight_upper + \
                Tprofiles[lower, range(Tprofiles.shape[1])] * \
                (1-weight_upper)
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
                weight_upper + \
                Sprofiles[lower, range(Sprofiles.shape[1])] * \
                (1-weight_upper)

        return T0, S0
    
    def vertical_buoyancy(self):
        """Move particles vertically according to their buoyancy"""
        in_ocean = np.where(self.elements.z<=0)[0]
        if len(in_ocean) > 0:
            self.elements.z[in_ocean] = np.minimum(0,
                self.elements.z[in_ocean] + self.elements.terminal_velocity[in_ocean] * self.time_step.total_seconds())

        # check for minimum height/maximum depth for each particle accouting also for
        # the sea surface height
        Zmin = -1.*(self.environment.sea_floor_depth_below_sea_level + self.environment.sea_surface_height)

        # Let particles stick to bottom
        bottom = np.where(self.elements.z < Zmin)
        if len(bottom[0]) > 0:
            logger.debug('%s elements reached seafloor, interacting with bottom' % len(bottom[0]))
            self.interact_with_seafloor()
            self.bottom_interaction(Zmin)
    
    def update(self):
        """Update positions and properties of elements."""
        self.light_along_trajectory()

        self.water_column_stretching()

        # Simply move particles with ambient current
        self.advect_ocean_current()

        # Advect particles due to surface wind drag,
        # according to element property wind_drift_factor
        self.advect_wind()

        # Stokes drift
        self.stokes_drift()

        # Turbulent Mixing
        self.update_terminal_velocity()
        self.detect_salinity()
        if self.get_config('drift:vertical_mixing') is True:
            self.vertical_mixing()
        else:
            self.vertical_buoyancy()
        # Vertical advection
        self.vertical_advection()
        
        if self.get_config('general:deac') is True:
            self.deac()

