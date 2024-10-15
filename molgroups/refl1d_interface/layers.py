from typing import List, Literal, Union, Optional
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import trapezoid

from refl1d.names import Parameter, SLD, Stack, Slab
from refl1d.model import Layer

from .groups import BaseGroupInterface, MolgroupsInterface

@dataclass(init=False)
class MolgroupsLayer(Layer):

    base_group: BaseGroupInterface
    normarea_group: MolgroupsInterface | None = None
    add_groups: List[MolgroupsInterface] = field(default_factory=list)
    overlay_groups: List[MolgroupsInterface] = field(default_factory=list)
    contrast: SLD | None=None
    thickness: float | Parameter = 0.0
    name: str | None = None

    def __init__(self,
                 base_group: BaseGroupInterface,
                 normarea_group: MolgroupsInterface | None = None,
                 add_groups: List[MolgroupsInterface] = [],
                 overlay_groups: List[MolgroupsInterface] = [],
                 contrast: SLD | None=None,
                 thickness: float | Parameter = 0.0,
                 name=None) -> None:

        if not isinstance(base_group, BaseGroupInterface):
            raise TypeError(f'Base group {base_group} must be an instance of BaseGroupInterface')

        self.base_group = base_group
        self.normarea_group = normarea_group
        self.add_groups = add_groups
        self.overlay_groups = overlay_groups

        if name is None:
            name = self.base_group.name

        self.thickness = Parameter.default(thickness, name=name+" thickness")
        self.interface = Parameter.default(0.0, name=name+" interface")

        self.name = name
        self.magnetism = None
        self.contrast = contrast

        self._penalty = 0.0

    def update(self):

        normarea: float | None = None

        # 0. update normarea_group

        # 1. update base_group (may be required twice if normarea_group is defined differently)
        self.base_group.update(bulknsld=self.contrast.rho.value)
        if self.normarea_group is None:
            normarea = self.base_group.normarea.value
        else:
            self.normarea_group.update(bulknsld=self.contrast.rho.value)
            if not hasattr(self.normarea_group, 'normarea'):
                print(f'Warning: {self.normarea_group.name} does not have normarea and cannot be normarea_group. Ignoring.')
            else:
                normarea = self.normarea_group.normarea.value
                if hasattr(self.base_group._molgroup, 'fnSetNormarea') & (self.base_group != self.normarea_group):
                    self.base_group._molgroup.fnSetNormarea(normarea)                
                self.base_group.normarea.value = normarea
                self.base_group.update(bulknsld=self.contrast.rho.value)

        # 2. apply normarea to all remaining objects
        for group in self.add_groups + self.overlay_groups:
            if hasattr(group._molgroup, 'fnSetNormarea') & (group != self.normarea_group):
                group._molgroup.fnSetNormarea(normarea)
        
        # 3. update all remaining objects
        for group in self.add_groups + self.overlay_groups:
            if group != self.normarea_group:
                group.update(bulknsld=self.contrast.rho.value)

    def profile(self, z):

        # 1. Write base group
        normarea, area, nsl = self.base_group.render(z)

        # 2. Add up all add_groups
        for group in self.add_groups:
            _, newarea, newnsl = group.render(z)
            area += newarea
            nsl += newnsl

        # 3. Add up all overlay_groups
        overlay_area = np.zeros_like(area)
        overlay_nsl = np.zeros_like(nsl)
        for group in self.overlay_groups:
            _, newarea, newnsl = group.render(z)
            overlay_area += newarea
            overlay_nsl += newnsl       

        # 4. Perform overlay
        # 4a. Calculate fraction of base_group + add_groups left after replacement
        frac_replacement = np.ones_like(area)
        if len(self.overlay_groups):
            over_filled = (area + overlay_area) > normarea
            frac_replacement[over_filled] = (area / (normarea - overlay_area))[over_filled]
        
        # 4b. Scale area, nsl by replacement fraction and add
        # in overlay_area, overlay_nsl
        area /= frac_replacement
        nsl /= frac_replacement
        area += overlay_area
        nsl += overlay_nsl

        # 4c. Scale stored profiles by frac_replacement
        for group in [self.base_group] + self.add_groups:
            group._stored_profile['frac_replacement'] = frac_replacement

        return normarea, area, nsl        

    def parameters(self):

        # TODO: figure out how to get all the unique parameters from the sub-objects
        # This will probably break if there are overlapping parameters

        return {k : p
                for group in [self.base_group] + self.add_groups + self.overlay_groups
                    for k, p in group._get_parameters().items()}


    def penalty(self) -> float:

        return self._penalty

    def _filled_profile(self, z):
        """Given area and nSL profiles, fill in the remaining volume with bulk material"""
        
        self.update()
        normarea, area, nsl = self.profile(z)

        # calculate penalty due to overfilling anywhere
        over_filled = area > normarea
        self._penalty = trapezoid((area - normarea)[over_filled], z[over_filled])

        # Fill in the remaining volume with buffer of appropriate nSLD
        nsld = 1e6 * nsl / (normarea * np.gradient(z)) + (1.0 - area / normarea) * self.contrast.rho.value

        # Return nSLD profile in Refl1D units
        return nsld

    def render(self, probe, slabs) -> None:
        """Adapted from refl1d.flayer.FunctionalProfile
        """
        Pw, Pz = slabs.microslabs(self.thickness)
        if len(Pw) == 0:
            return

        # add molgroups layer slabs
        slabs.extend(rho=[self._filled_profile(Pz)], irho=[np.zeros_like(Pz)], w=Pw)

        # add buffer slab (better done at sample stack level)
        #slabs.append(rho=self.contrast.rho.value, irho=0.0, w=0.0, sigma=0.0)

    def _get_moldat(self) -> dict:
        """Gets moldat object for plotting and statistical analysis

        Returns:
            dict: moldat
        """

        return {group.name: group._stored_profile
                    for group in [self.base_group] + self.add_groups + self.overlay_groups}

# =============

@dataclass(init=False, eq=False, match_args=False)
class MolgroupsStack(Stack):

    substrate: Stack | Slab
    molgroups_layer: MolgroupsLayer

    def __init__(self,
                 substrate: Stack | Slab,
                 molgroups_layer: MolgroupsLayer,
                 name="MolgroupsStack",
                 **kw):
        layer_contrast = Slab(material=molgroups_layer.contrast, thickness=0.0000, interface=0.0000)
        #if isinstance(substrate, Stack):
        #    substrate.layers[-1].thickness = substrate.layers[-1].thickness - molgroups_layer.base_group.overlap
        #elif isinstance(substrate, Slab):
        #    substrate.thickness = substrate.thickness - molgroups_layer.overlap
        layers = substrate | molgroups_layer | layer_contrast
        super().__init__(base=None, 
                         layers=layers,
                         name=name,
                         interface=None,
                         thickness=None)
        
        self.substrate, self.molgroups_layer = substrate, molgroups_layer

    def __copy__(self):
        stack = MolgroupsStack(substrate=self.substrate,
                               molgroups_layer=self.molgroups_layer)
        stack.interface = self.interface
        stack._layers = self._layers[:]
        stack.thickness = self.thickness

        return stack

    def __getstate__(self):
        return self.interface, self._layers, self.name, self.thickness, self.substrate, self.molgroups_layer

    def __setstate__(self, state):
        self.interface, self._layers, self.name, self.thickness, self.substrate, self.molgroups_layer = state
        # TODO: not clear that this is needed here.  The thickness parameter
        # from __getstate__ should have a valid expression in it.
        self._set_thickness()