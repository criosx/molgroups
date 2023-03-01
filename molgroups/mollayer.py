import numpy

from refl1d.model import Layer
from bumps.parameter import Parameter, flatten
from refl1d.util import merge_ends

from molgroups.mol import CompositenSLDObj


class MolLayer(Layer):
    """Subclass of refl1d Layer base object for use with molgroups layers"""
    def __init__(self, bulknsld, parameter_map, base_groups=[], overlay_groups=[], thickness=None,
                 interface=None, tol=1e-3, name=None) -> None:
        self.name = name
        self.thickness = Parameter.default(thickness, name=name+" thickness")
        self.interface = Parameter.default(interface, name=name+" interface")
        self.base_groups = base_groups
        self.overlay_groups=overlay_groups
        self.bulknsld = bulknsld
        self.tol = tol
        self.parameter_map = parameter_map

        if not len(base_groups):
            raise ValueError(f'At least one base group required in {name}')
        elif len(base_groups) == 1:
            self._base = base_groups[0]
        else:
            self._base = CompositenSLDObj(base_groups, name=f'{name}.base')
        
        if not len(overlay_groups):
            self._overlay = None
        elif len(overlay_groups) == 1:
            self._overlay = overlay_groups[0]
        else:
            self._overlay = CompositenSLDObj(overlay_groups, name=f'{name}.overlay')

    def parameters(self):
        # Find all refl1d/bumps parameters and return dict
        #pars = self._get_parameters()

        #return {par.name: par for par in pars}

        pars = flatten([p.parameters() if not isinstance(p, list) else [item.parameters() for item in p] for _, _, p in self.parameter_map])

        return {p.name: p for p in pars}

    def update_molgroups(self):

        for gp, attr, obj in self.parameter_map:
            if isinstance(obj, list):
                setattr(gp, attr, [float(item) for item in obj])
            else:
                setattr(gp, attr, obj.value)

        self._base.update()

    def profile(self, z):
        """Calculates total scattering length density profile from bulk nSLD"""

        # calculate total area from base groups
        area, nsl, _ = self._base.fnGetProfiles(z)
        normarea = area.max()

        # calculate total area from proteins and overlay proteins on bilayers
        if self._overlay is not None:
            area, nsl = self._overlay.fnOverlayProfile(z, area, nsl, normarea)

        # Fill in the remaining volume with buffer of appropriate nSLD
        nsld = nsl / (normarea * numpy.gradient(z)) + (1.0 - area / normarea) * float(self.bulknsld)

        # Return nSLD profile in Refl1D units
        return nsld * 1e6

    def render(self, probe, slabs):
        # refl1d renderer
        Pw, Pz = slabs.microslabs(self.thickness.value)
        if len(Pw) == 0:
            return

        # Calculate 
        self.update_molgroups()
        phi = numpy.asarray(self.profile(Pz))
        if phi.shape != Pz.shape:
            raise TypeError("profile function '%s' did not return array phi(z)"
                            %self.__name__)
        Pw, phi = merge_ends(Pw, phi, tol=self.tol)
        #P = M*phi + S*(1-phi)
        slabs.extend(rho=[numpy.real(phi)], irho=[numpy.imag(phi)], w=Pw)
