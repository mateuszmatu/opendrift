"""    def premature_deactivation(self):

        indices = []        
        if self.get_config('general:premature_deactivation') is None:
            logger.info('Deactivation function deactivated')
            return

        # Hard limits for temperature or salinity (basic usage)
        elif self.get_config('general:premature_deactivation') is 'temperature_minmax':
            indices = [el < self.get_config('general:deactivation_min') or el > self.get_config('general:deactivation_max') for el in self.environment.sea_water_temperature]

        elif self.get_config('general:premature_deactivation') is 'salinity_minmax':
            indices = [el < self.get_config('general:deactivation_min') or el > self.get_config('general:deactivation_max') for el in self.environment.sea_water_salinity]

        elif self.get_config('general:premature_deactivation') is 'shortwave_minmax':
            Qdown = self.shortwave_radiation_at_depth(self.environment.net_downward_shortwave_flux_at_sea_water_surface, self.elements.z)
            indices = [el < self.get_config('general:deactivation_min') or el > self.get_config('general:deactivation_max') for el in Qdown]
            
        # Uses the new health attribute of pelagic egg.
        # Eg. 1 particle is an egg population, how much % of population survives along trajectory
        # Eg. 2 Likelyhood of egg surviving trajectory
        elif self.get_config('general:premature_deactivation') is 'exposure':
            # Create some exposure "health drain" functions
            if self.get_config('general:deactivation_exposure') is 'simple_minmax_exposure':
                health_indices = [el < self.get_config('general:deactivation_min') or 
                           el > self.get_config('general:deactivation_max') 
                           for el in self.environment.sea_water_temperature]

                self.elements.health[health_indices == np.True_] -= 20
            
            elif self.get_config('general:deactivation_exposure') is 'randomly_consumed':
                self.elements.health -= np.random.rand(len(self.elements.health)) * 10

            indices = [el <= 0 for el in self.elements.health]

        if len(indices) > 0:
            self.deactivate_elements(indices, 'Deactivated through general:premature_deactivation config.')

                'The fetch length when using tabularised Stokes drift.'},
'general:seafloor_action': {'type': 'enum', 'default': 'lift_to_seafloor',
    'enum': ['none', 'lift_to_seafloor', 'deactivate', 'previous'],
    'description': '"deactivate": elements are deactivated; "lift_to_seafloor": elements are lifted to seafloor level; "previous": elements are moved back to previous position; "none"; seafloor is ignored.',
    'level': CONFIG_LEVEL_ADVANCED},
'drift:truncate_ocean_model_below_m': {'type': 'float', 'default': None,
    'min': 0, 'max': 10000, 'units': 'm',
    'description': 'Ocean model data are only read down to at most this depth, and extrapolated below. May be specified to read less data to improve performance.',
    'level': CONFIG_LEVEL_ADVANCED},
    'seed:z': {'type': 'float', 'default': 0,
        'min': -10000, 'max': 0, 'units': 'm',
    'description': 'Depth below sea level where elements are released. This depth is neglected if seafloor seeding is set selected.',
    'level': CONFIG_LEVEL_ESSENTIAL},
'seed:seafloor': {'type': 'bool', 'default': False,
    'description': 'Elements are seeded at seafloor, and seeding depth (z) is neglected.',
    'level': CONFIG_LEVEL_ESSENTIAL},
'general:premature_deactivation': {'type': 'enum', 'default': None, 'enum': [None, 'temperature_minmax', 'salinity_minmax', 'exposure', 'shortwave_minmax'],
    'description': 'Deactivation function', 
    'level': CONFIG_LEVEL_ADVANCED},
'general:deactivation_min': {'type': 'float', 'default': -999.0, 'min': -999, 'max': 999, 'units': 'None',
    'description':  'Min thresholds for the deactivation function.',
    'level': CONFIG_LEVEL_ADVANCED},
'general:deactivation_max': {'type': 'float', 'default': 999.0, 'min': -999, 'max': 999, 'units': 'None',
    'description':  'Max thresholds for the deactivation function.',
    'level': CONFIG_LEVEL_ADVANCED},
'general:deactivation_exposure': {'type': 'enum', 'default': 'simple_minmax_exposure', 
    'enum': ['simple_minmax_exposure', 'randomly_consumed'],
    'description': 'Exposure functions for premature_deactivation.',
    'level': CONFIG_LEVEL_ADVANCED}
})

from opendrift.models.pelagicegg import PelagicEggDrift, Lagrangian3DArray
import numpy as np


class PelagicEggHealth(Lagrangian3DArray):


    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # for NEA Cod
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': '[]',
                                       'default': 31.25}),  # for NEA Cod
        ('density', {'dtype': np.float32,
                     'units': 'kg/m^3',
                     'default': 1028.}),
        ('hatched', {'dtype': np.float32,
                     'units': '',
                     'default': 0.}),
        ('health', {'dtype': np.float32,
                     'units': '',
                     'default': 100.})])


class PelagicEggHealthDrift(PelagicEggDrift):

    ElementType = PelagicEggHealth
    required_variables = {
    'x_sea_water_velocity': {'fallback': 0},
    'y_sea_water_velocity': {'fallback': 0},
    'sea_surface_height': {'fallback': 0},
    'sea_surface_wave_significant_height': {'fallback': 0},
    'sea_ice_area_fraction': {'fallback': 0},
    'x_wind': {'fallback': 0},
    'y_wind': {'fallback': 0},
    'land_binary_mask': {'fallback': None},
    'sea_floor_depth_below_sea_level': {'fallback': 100},
    'ocean_vertical_diffusivity': {'fallback': 0.02, 'profiles': True},
    'ocean_mixed_layer_thickness': {'fallback': 50},
    'sea_water_temperature': {'fallback': 10, 'profiles': True},
    'sea_water_salinity': {'fallback': 34, 'profiles': True},
    'surface_downward_x_stress': {'fallback': 0},
    'surface_downward_y_stress': {'fallback': 0},
    'turbulent_kinetic_energy': {'fallback': 0},
    'turbulent_generic_length_scale': {'fallback': 0},
    'upward_sea_water_velocity': {'fallback': 0},
    'net_downward_shortwave_flux_at_sea_water_surface': {'fallback': 145}
    }


    # Default colors for plotting
    status_colors = {'initial': 'green', 'active': 'blue',
                     'hatched': 'red', 'eaten': 'yellow', 'died': 'magenta'}

    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(PelagicEggDrift, self).__init__(*args, **kwargs)"""


import numpy as np
import logging; logger = logging.getLogger(__name__)

from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

class MyElement(Lagrangian3DArray):
    """Extending Lagrangian3DArray with specific properties for pelagic eggs
    """

    variables = Lagrangian3DArray.add_variables([
        ('health_percentage', {'dtype': np.float32,
                     'units': '',
                     'default': 100.}),
        ('light', {'dtype': np.float32,
                   'units': 'W m^-2',
                   'default': 0})
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
        'net_downward_shortwave_flux_at_sea_water_surface': {'fallback': 0}
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
            'general:deac': {'type': 'enum', 'default': None, 'enum': [None, 'temperature_minmax', 'salinity_minmax', 'exposure', 'shortwave_minmax'],
                'description': 'Deactivation function', 
                'level': CONFIG_LEVEL_ADVANCED},
            'general:deac_min': {'type': 'float', 'default': -999.0, 'min': -999, 'max': 999, 'units': 'None',
                'description':  'Min thresholds for the deactivation function.',
                'level': CONFIG_LEVEL_ADVANCED},
            'general:deac_max': {'type': 'float', 'default': 999.0, 'min': -999, 'max': 999, 'units': 'None',
                'description':  'Max thresholds for the deactivation function.',
                'level': CONFIG_LEVEL_ADVANCED},
            'general:deac_exposure': {'type': 'enum', 'default': 'simple_minmax_exposure', 
                'enum': ['simple_minmax_exposure', 'randomly_consumed'],
                'description': 'Exposure functions for premature_deactivation.',
                'level': CONFIG_LEVEL_ADVANCED}
            })

    def deac(self):

        indices = []        
        if self.get_config('general:deac') is None:
            logger.info('Deactivation function deactivated')
            return

        # Hard limits for temperature or salinity (basic usage)
        elif self.get_config('general:deac') is 'temperature_minmax':
            indices = [el < self.get_config('general:deac_min') or el > self.get_config('general:deac_max') for el in self.environment.sea_water_temperature]

        elif self.get_config('general:deac') is 'salinity_minmax':
            indices = [el < self.get_config('general:deac_min') or el > self.get_config('general:deac_max') for el in self.environment.sea_water_salinity]

        elif self.get_config('general:deac') is 'shortwave_minmax':
            Qdown = self.shortwave_radiation_at_depth(self.environment.net_downward_shortwave_flux_at_sea_water_surface, self.elements.z)
            print(Qdown)
            indices = [el < self.get_config('general:deac_min') or el > self.get_config('general:deac_max') for el in Qdown]
            
        # Uses the new health attribute of pelagic egg.
        # Eg. 1 particle is an egg population, how much % of population survives along trajectory
        # Eg. 2 Likelyhood of egg surviving trajectory
        elif self.get_config('general:deac') is 'exposure':
            # Create some exposure "health drain" functions
            if self.get_config('general:deac_exposure') is 'simple_minmax_exposure':
                health_indices = [el < self.get_config('general:deac_min') or 
                           el > self.get_config('general:deac_max') 
                           for el in self.environment.sea_water_temperature]

                self.elements.health[health_indices == np.True_] -= 20
            
            elif self.get_config('general:deac_exposure') is 'randomly_consumed':
                self.elements.health -= np.random.rand(len(self.elements.health)) * 10

            indices = [el <= 0 for el in self.elements.health]

        if len(indices) > 0:
            self.deactivate_elements(indices, 'Deactivated through general:premature_deactivation config.')

    def update(self):
        """Update positions and properties of elements."""

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
        else:  # Buoyancy
            self.vertical_buoyancy()

        # Vertical advection
        self.vertical_advection()
        print('running deac')
        self.deac()

