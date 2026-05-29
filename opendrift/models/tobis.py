# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2021, Trond Kristiansen, Niva
# Jan 2021 Simplified by Knut-Frode Dagestad, MET Norway, and adapted to to Kvile et al. (2018)

import datetime
import numpy as np
from scipy.interpolate import RectBivariateSpline
import logging; logger = logging.getLogger(__name__)
from opendrift.models.larvalfish import LarvalFishElement
from opendrift.models.oceandrift import Lagrangian3DArray, OceanDrift
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED


class TobisElement(Lagrangian3DArray):
    """
    Extending Lagrangian3DArray with specific properties for Tobis
    """

    variables = Lagrangian3DArray.add_variables([
        ('diameter', {'dtype': np.float32,
                      'units': 'm',
                      'default': 0.0014}),  # need adjustment for Tobis if egg advection is turned on
        ('neutral_buoyancy_salinity', {'dtype': np.float32,
                                       'units': 'PSU',
                                       'default': 31.25}),  # need adjustment for Tobis if egg advection is turned on
        ('stage_fraction', {'dtype': np.float32,  # to track percentage of development time completed
                            'units': '',
                            'default': 0.}),
        ('phase', {'dtype': np.uint8,  # 0 for eggs, 1 for larvae, 2 for fish
                     'units': '',
                     'default': 0}),
        ('hatch_rate', {'dtype': np.float32,
                    'units': '',
                    'default': 0}),
        ('length', {'dtype': np.float32,
                    'units': 'mm',
                    'default': 0}),
        ('weight', {'dtype': np.float32,
                    'units': 'mg',
                    'default': 0.08})])


class Tobis(OceanDrift):
    """Buoyant particle trajectory model based on the OpenDrift framework.

        Developed at MET Norway

        Generic module for particles that are subject to vertical turbulent
        mixing with the possibility for positive or negative buoyancy

        Particles could be e.g. oil droplets, plankton, or sediments

    """

    ElementType = TobisElement

    required_variables = {
        'x_sea_water_velocity': {'fallback': 0},
        'y_sea_water_velocity': {'fallback': 0},
        'sea_surface_height': {'fallback': 0},
        'sea_surface_wave_significant_height': {'fallback': 0},
        'x_wind': {'fallback': 0},
        'y_wind': {'fallback': 0},
        'land_binary_mask': {'fallback': None},
        'sea_floor_depth_below_sea_level': {'fallback': 100},
        'ocean_vertical_diffusivity': {'fallback': 0.01, 'profiles': True},
        'ocean_mixed_layer_thickness': {'fallback': 50},
        'sea_water_temperature': {'fallback': 10, 'profiles': True},
        'sea_water_salinity': {'fallback': 34, 'profiles': True},
        'sea_surface_wave_stokes_drift_x_velocity': {'fallback': 0},
        'sea_surface_wave_stokes_drift_y_velocity': {'fallback': 0},
    }


    def __init__(self, *args, **kwargs):

        # Calling general constructor of parent class
        super(Tobis, self).__init__(*args, **kwargs)

        # IBM configuration options
        self._add_config({
            'IBM:fraction_of_timestep_swimming':
                {'type': 'float', 'default': 0.15,
                 'min': 0.0, 'max': 1.0, 'units': 'fraction',
                 'description': 'Fraction of timestep swimming',
                 'level': CONFIG_LEVEL_ADVANCED},
            })
        
        self._add_config({
            'drift:egg_advection':
            {'type': 'bool', 'default':False,
             'description':'Turn on/off egg advection',
             'level':CONFIG_LEVEL_ADVANCED}
        })

        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:vertical_mixing_at_surface', True)
        self._set_config_default('drift:vertical_advection_at_surface', True)

        self.hatch_time = self.get_hatch_time_func()

    def get_hatch_time_func(self):
        """
        Total hatch time based on Smigielski et al. (1984), doi: 10.3354/meps014287

        The underlying data table is taken directly from the paper. Interpolation is
        linear with date and second order in the rate direction.
        """

        days_tab = np.array([
            [61, 51, 39, 25],
            [82, 67, 48, 30],
            [135, 116, 82, 55],
        ])
        temp_tab = np.array([2, 4, 7, 10])
        rate_tab = np.array([0, 0.5, 1])

        spline = RectBivariateSpline(
            x=rate_tab,
            y=temp_tab,
            z=days_tab,
            kx=2,
            ky=1,
        )

        def hatch_time_fn(rate, temp):
            """
            Total hatch time based on Smigielski et al. (1984), doi: 10.3354/meps014287

            The underlying data table is taken directly from the paper. Interpolation is
            linear with date and second order in the rate direction.

            :param rate: 0 = earliest spawners, 0.5 = median spawners, 1 = latest spawners
            :param temp: Ambient temperature, in degrees Celcius
            :return: Total hatch time, in days
            """
            temp = np.minimum(temp_tab[-1], np.maximum(temp_tab[0], temp))
            out = spline(rate.ravel(), temp.ravel(), grid=False)
            return out.reshape(rate.shape)

        return hatch_time_fn

    def update_terminal_velocity(self, Tprofiles=None,
                                 Sprofiles=None, z_index=None):
        """Calculate terminal velocity for Pelagic Egg

        according to
        S. Sundby (1983): A one-dimensional model for the vertical
        distribution of pelagic fish eggs in the mixed layer
        Deep Sea Research (30) pp. 645-661

        Method copied from ibm.f90 module of LADIM:
        Vikebo, F., S. Sundby, B. Aadlandsvik and O. Otteraa (2007),
        Fish. Oceanogr. (16) pp. 216-228
        """
        g = 9.81  # ms-2

        # Pelagic Egg properties that determine buoyancy
        eggsize = self.elements.diameter  # 0.0014 for NEA Cod
        eggsalinity = self.elements.neutral_buoyancy_salinity
        # 31.25 for NEA Cod

        # prepare interpolation of temp, salt
        if not (Tprofiles is None and Sprofiles is None):
            if z_index is None:
                z_i = range(Tprofiles.shape[0])  # evtl. move out of loop
                # evtl. move out of loop
                z_index = interp1d(-self.environment_profiles['z'],
                                   z_i, bounds_error=False)
            zi = z_index(-self.elements.z)
            upper = np.maximum(np.floor(zi).astype(np.uint8), 0)
            lower = np.minimum(upper + 1, Tprofiles.shape[0] - 1)
            weight_upper = 1 - (zi - upper)

        # do interpolation of temp, salt if profiles were passed into
        # this function, if not, use reader by calling self.environment
        if Tprofiles is None:
            T0 = self.environment.sea_water_temperature
        else:
            T0 = Tprofiles[upper, range(Tprofiles.shape[1])] * \
                 weight_upper + \
                 Tprofiles[lower, range(Tprofiles.shape[1])] * \
                 (1 - weight_upper)
        if Sprofiles is None:
            S0 = self.environment.sea_water_salinity
        else:
            S0 = Sprofiles[upper, range(Sprofiles.shape[1])] * \
                 weight_upper + \
                 Sprofiles[lower, range(Sprofiles.shape[1])] * \
                 (1 - weight_upper)

        # The density difference between a pelagic egg and the ambient water
        # is regulated by their salinity difference through the
        # equation of state for sea water.
        # The Egg has the same temperature as the ambient water and its
        # salinity is regulated by osmosis through the egg shell.
        DENSw = self.sea_water_density(T=T0, S=S0)
        DENSegg = self.sea_water_density(T=T0, S=eggsalinity)
        dr = DENSw - DENSegg  # density difference

        # water viscosity
        my_w = 0.001 * (1.7915 - 0.0538 * T0 + 0.007 * (T0 ** (2.0)) - 0.0023 * S0)
        # ~0.0014 kg m-1 s-1

        # terminal velocity for low Reynolds numbers
        W = (1.0 / my_w) * (1.0 / 18.0) * g * eggsize ** 2 * dr

        # check if we are in a Reynolds regime where Re > 0.5
        highRe = np.where(W * 1000 * eggsize / my_w > 0.5)

        # Use empirical equations for terminal velocity in
        # high Reynolds numbers.
        # Empirical equations have length units in cm!
        my_w = 0.01854 * np.exp(-0.02783 * T0)  # in cm2/s
        d0 = (eggsize * 100) - 0.4 * \
             (9.0 * my_w ** 2 / (100 * g) * DENSw / dr) ** (1.0 / 3.0)  # cm
        W2 = 19.0 * d0 * (0.001 * dr) ** (2.0 / 3.0) * (my_w * 0.001 * DENSw) ** (-1.0 / 3.0)
        # cm/s
        W2 = W2 / 100.  # back to m/s

        W[highRe] = W2[highRe]
        self.elements.terminal_velocity = W
    
    
    def update_egg_stage_fraction(self):
        """Egg development according to Christensen et al. (2008), doi:10.1139/F08-073

        Using initialization model `e` (variable maturation)"""

    
        eggs = np.where(self.elements.phase==0)[0]
        if len(eggs) > 0:
            dev_days = self.hatch_time(self.elements.hatch_rate[eggs], self.environment.sea_water_temperature[eggs])
            stage_increase = self.time_step.total_seconds()/(dev_days*24*60*60)
            self.elements.stage_fraction[eggs] += stage_increase
            hatching = np.where(self.elements.stage_fraction[eggs]>=1)[0]
            if len(hatching) > 0:
                logger.debug('Hatching %s eggs' % len(hatching))
                self.elements.phase[eggs[hatching]] = 1  # Eggs with stage_fraction >= 1 are hatched (1)
        
        
    def update_larval_stage_fraction(self):
        """Larval development according to Christensen et al. (2008), doi:10.1139/F08-073

        The function modifies (in-place) the variable `stage`. The length of the larvae is
        assumed to be related to the stage as

        L = L0 * (1 - s) + Lm * s

        where L0 = 7.73mm, Lm = 40mm and s = stage - 1"""
        
        larvae = np.where(self.elements.phase==1)[0]

        if len(larvae) > 0:
            s = self.elements.stage_fraction[larvae] - 1
            Lm = 40
            L0 = 7.73
            L_inf = 218
            
            L = L0 + s * (Lm - L0)
            gam = 0.316
            lamb0 = -1.725
            lamb1 = 0.136
            lamb = np.exp(lamb0 + lamb1 * self.environment.sea_water_temperature[larvae])

            dLdt = lamb * np.power(L / L0, gam) * (1 - L / L_inf)
            L_new = L + dLdt * self.time_step.total_seconds()/86400
            self.elements.stage_fraction[larvae] = 1 + (L_new - L0) / (Lm - L0)
            self.elements.length[larvae] = L_new

            metamorphing = np.where(self.elements.stage_fraction[larvae] >= 2)[0]
            if len(metamorphing) > 0:
                logger.info('Metamorphing %s larvae' % len(metamorphing))
                self.elements.phase[larvae[metamorphing]] = 2  # Larvae with stage_fraction >= 2 are metamorphing (2)
                self.deactivate_metamorphed()

        ### TODO: missing a method for increasing weight. Could use whatever was in OpenDrift. 

    def larvae_vertical_migration(self):

        larvae = np.where(self.elements.phase==1)[0]
        if len(larvae) == 0:
            return

        # Vertical migration of Larvae
        # Swim function from Peck et al. 2006
        L = self.elements.length[larvae]
        swim_speed = (0.261*(L**(1.552*L**(-0.08))) - 5.289/L) / 1000  # TODO: this was here previously
        f = self.get_config('IBM:fraction_of_timestep_swimming')
        max_migration_per_timestep = f*swim_speed*self.time_step.total_seconds()

        # Using here UTC hours. Should be changed to local solar time,
        # although a phase shift of some hours should not make much difference
        if self.time.hour < 12:
            direction = -1  # Swimming down when light is increasing
        else:
            direction = 1  # Swimming up when light is decreasing

        self.elements.z[larvae] = np.minimum(0, self.elements.z[larvae] + direction*max_migration_per_timestep)
    
    def freeze_egg(self):
        eggs = np.where(self.elements.phase==0)[0]
        larvae = np.where(self.elements.phase==1)[0]
        if len(eggs) > 0:
            self.elements.moving[eggs] = 0
        if len(larvae) > 0:
            self.elements.moving[larvae] = 1
        
    def deactivate_metamorphed(self):
        metamorphed_indices = [el >= 2 for el in self.elements.phase]            
        self.deactivate_elements(metamorphed_indices, reason='Metamorphosed')

    def update(self):

        self.update_egg_stage_fraction()
        self.update_larval_stage_fraction()
        if self.get_config('drift:egg_advection') is False:
            self.freeze_egg()
        self.advect_ocean_current()

        # Stokes drift
        self.stokes_drift()

        self.update_terminal_velocity()
        self.vertical_mixing()
        self.larvae_vertical_migration()
