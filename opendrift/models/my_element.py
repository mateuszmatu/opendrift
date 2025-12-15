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
                           'units': 'W m^-2',
                           'default': 0}),
        ('vertical_swim_speed', {'dtype': np.float32,
                                  'units': 'm s^-1',
                                  'default': 0}),
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # for NEA Cod
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': '[]',
                                       'default': 31.25}),  # for NEA Cod
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
        'net_downward_shortwave_flux_at_sea_water_surface': {'fallback': 100}
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
                                  'level': CONFIG_LEVEL_ADVANCED}
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
        
        if self.elements.light > self.elements.prefered_light:
            W = -self.elements.vertical_swim_speed
        elif self.elements.light < self.elements.prefered_light:
            W = self.elements.vertical_swim_speed

        self.elements.terminal_velocity = W


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
        if self.get_config('drift:vertical_mixing') is True:
            self.vertical_mixing()
        else:
            self.vertical_buoyancy()
        # Vertical advection
        self.vertical_advection()
        
        if self.get_config('general:deac') is True:
            self.deac()

