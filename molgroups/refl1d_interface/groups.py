"""Refl1D shims for molgroups.mol objects.
    Uses Refl1D Parameter objects for parameters, thus allowing model serialization
"""

from typing import List, Tuple, Callable, Dict, TypedDict
from dataclasses import dataclass, field, fields
import uuid
import functools

import numpy as np

from scipy.integrate import trapezoid
from refl1d.names import Parameter
from bumps.parameter import Calculation

from molgroups.mol import nSLDObj, BLM, ssBLM, tBLM, Box2Err, BoxHermite, BLMProteinComplex
from molgroups.components import Component, Lipid, Tether, bme

from periodictable.fasta import H2O_SLD, D2O_SLD

def sld_from_bulk(rhoH: float, rhoD: float, bulknsld: float) -> float:
    """Calculates scattering length density of material with
        labile hydrogens from bulk nSLD and SLDs in pure water and D2O

    Args:
        rhoH (float): nSLD in pure H2O
        rhoD (float): nSLD in pure D2O
        bulknsld (float): nSLD of bulk water

    Returns:
        float: nSLD of material 
    """
    
    frac_d2o = (bulknsld - H2O_SLD) / (D2O_SLD - H2O_SLD)

    return rhoH + (rhoD - rhoH) * frac_d2o

class ReferencePoint(Parameter):

    def __init__(self, function: Callable | None = None, description: str = '', name: str | None = None, id: str | None = None, discrete: bool = False, tags: List[str] | None = None, **kw):
        calculation = Calculation(description=description)
        if function is not None:
            calculation.set_function(function)
        tags = [] if tags is None else tags
        super().__init__(slot=calculation, fixed=True, name=name, id=id, discrete=discrete, tags=tags + ['Reference Point'], **kw)
    
    def set_function(self, function: Callable) -> None:

        self.slot.set_function(function)

@dataclass
class MolgroupsInterface:
    """Base class for interacting with molgroups objects
    """

    id: str | None = None
    name: str | None = None
    nf: Parameter = field(default_factory=lambda: Parameter(name='number fraction', value=1))
    _molgroup: nSLDObj | None = None
    _stored_profile: dict | None = None
    _group_names: dict[str, List[str]] = field(default_factory=dict)

    def __post_init__(self) -> None:

        if self.id is None:
            self.id = str(uuid.uuid4())

        if not self._group_names:
            self._group_names = {f'{self.name}': [f'{self.name}']}

        for f in fields(self):
            if f.type == Parameter:
                default_name = f.default_factory().name
                p = getattr(self, f.name)
                if hasattr(p, 'name'):
                    if p.name == default_name:
                        p.name = f'{self.name} {p.name}'
                        setattr(self, f.name, p)
                else:
                    setattr(self, f.name, Parameter.default(p, name=f'{self.name} {default_name}'))

    def _get_parameters(self) -> dict[str, Parameter]:
        """Gets a list of the parameters associated with the interactor

        Returns:
            List[Parameter]: Parameter list
        """

        pars = {}
        for f in fields(self):
            if f.type == Parameter:
                p = getattr(self, f.name)
                p = Parameter.default(p, name=f'{self.name} {f.name}')
                setattr(self, f.name, p)
                pars.update({f'{self.name} {f.name}': p})
            elif f.type == List[Parameter]:
                plist = getattr(self, f.name)
                for i, p in enumerate(plist):
                    p = Parameter.default(p, name=f'{self.name} {f.name}{i}')
                    pars.update({f'{self.name} {f.name}{i}': p})
                    plist[i] = p
                setattr(self, f.name, plist)
            elif f.type == ReferencePoint:
                p: ReferencePoint = getattr(self, f.name)
                pars.update({f'{self.name} {p.name}': p})

        return pars

    def update(self, bulknsld: float) -> None:
        """Updates the molecular group with current values of the parameters,
            usually by calling fnSet
        """

        pass

    def old_render(self, z: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        """Renders the molecular group to an area and nSL

        Args:
            z (np.ndarray): spatial domain on which to render molecular group

        Returns:
            Tuple (float,np.ndarray, np.ndarray): normarea, area, nSL
        """

        normarea, area, nsl = self._molgroup.fnWriteProfile(z)

        return normarea, area, nsl
    
    def render(self, z: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        """Renders the molecular group to an area and nSL and stores result

        Args:
            z (np.ndarray): spatial domain on which to render molecular group

        Returns:
            Tuple (float,np.ndarray, np.ndarray): normarea, area, nSL
        """

        self.store_profile(z)
        area = self._stored_profile['area']
        nsl = self._stored_profile['sl']
        normarea = self._stored_profile['normarea']

        return normarea, area, nsl
    
    def store_profile(self, z: np.ndarray) -> dict:
        """Renders the molecular group and writes to a dict

        Args:
            z (np.ndarray): spatial domain on which to render molecular group

        Returns:
            dict: stored profile dictionary
        """

        self._stored_profile = self._molgroup.fnWriteGroup2Dict(dict(frac_replacement=1), self.name, z)
        self._stored_profile = self._molgroup.fnWriteProfile2Dict(self._stored_profile, z)
        self._stored_profile['normarea'] = self._stored_profile['area'].max()
    
    def set_group_names(self, group_names: dict[str, List[str]]) -> None:
        """Sets default group names

        Args:
            dict[str, List[str]]: group names structure
        """

        self._group_names = group_names

    def _center_of_volume(self):

        if self._stored_profile is None:
            return 0.0

        z, area = self._stored_profile['zaxis'], self._stored_profile['area']

        return trapezoid(area * z, z) / trapezoid(area, z) if np.sum(area) else 0.0

@dataclass
class BLMInterface(MolgroupsInterface):
    """Refl1D interactor for free floating BLM
    """

    _molgroup: BLM | None = None

    lipids: List[Lipid] = field(default_factory=list)
    inner_lipid_nf: List[Parameter] = field(default_factory=list)
    outer_lipid_nf: List[Parameter] = field(default_factory=list)
    startz: Parameter = field(default_factory=lambda: Parameter(name='position of inner hydrophobic interface', value=0.9))
    vf_bilayer: Parameter = field(default_factory=lambda: Parameter(name='volume fraction bilayer', value=0.9))
    l_lipid1: Parameter = field(default_factory=lambda: Parameter(name='inner acyl chain thickness', value=10.0))
    l_lipid2: Parameter = field(default_factory=lambda: Parameter(name='outer acyl chain thickness', value=10.0))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='bilayer roughness', value=5))
    normarea: Parameter = field(default_factory=lambda: Parameter(name='normarea', value=1))

    bilayer_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='bilayer_center', description='center of bilayer'))
    inner_headgroup_bottom: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_headgroup_bottom', description='bottom of inner headgroups'))
    inner_headgroup_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_headgroup_center', description='center of inner headgroups'))
    inner_hydrophobic_interface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_hydrophobic_interface', description='interface between inner headgroups and acyl chains'))
    outer_hydrophobic_interface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_hydrophobic_interface', description='interface between outer headgroups and acyl chains'))
    outer_headgroup_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_headgroup_center', description='center of outer headgroups'))
    outer_headgroup_top: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_headgroup_top', description='top of outer headgroups'))

    def __post_init__(self):
        self._molgroup = BLM(inner_lipids=self.lipids,
                             outer_lipids=self.lipids,
                             inner_lipid_nf=[p.value if hasattr(p, 'value') else p for p in self.inner_lipid_nf],
                             outer_lipid_nf=[p.value if hasattr(p, 'value') else p for p in self.outer_lipid_nf],
                             name=self.name)

        n_lipids = len(self.lipids)
        self._group_names = {f'{self.name} inner headgroups': [f'{self.name}.headgroup1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} inner acyl chains': [f'{self.name}.methylene1_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer acyl chains': [f'{self.name}.methylene2_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl2_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer headgroups': [f'{self.name}.headgroup2_{i}' for i in range(1, n_lipids + 1)],
                }

        # connect reference points
        self.bilayer_center.set_function(self._molgroup.fnGetCenter)
        self.inner_headgroup_bottom.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc - blm.av_hg1_l, self._molgroup))
        self.inner_headgroup_center.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc - 0.5 * blm.av_hg1_l, self._molgroup))
        self.inner_hydrophobic_interface.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc, self._molgroup))
        self.outer_hydrophobic_interface.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc, self._molgroup))
        self.outer_headgroup_center.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc + 0.5 * blm.av_hg2_l, self._molgroup))
        self.outer_headgroup_top.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc + blm.av_hg2_l, self._molgroup))

        super().__post_init__()

    def update(self, bulknsld: float):

        self._molgroup.fnSet(sigma=self.sigma.value,
            bulknsld=bulknsld * 1e-6,
            startz=self.startz.value,
            l_lipid1=self.l_lipid1.value,
            l_lipid2=self.l_lipid2.value,
            vf_bilayer=self.vf_bilayer.value,
            nf_inner_lipids=[p.value for p in self.inner_lipid_nf],
            nf_outer_lipids=[p.value for p in self.outer_lipid_nf],
            radius_defect=1e8)
        
        self.normarea.value = self._molgroup.normarea

# ============= BaseGroup objects ===============

@dataclass
class BaseGroupInterface(MolgroupsInterface):
    """Interface specifically for base groups, i.e. those that occupy the edges of the molgroups canvas
    """

    normarea: Parameter | float = 1.0
    overlap: Parameter | float = 20.0

    def __post_init__(self) -> None:

        self.normarea = Parameter.default(self.normarea, name=f'{self.name} normarea', fixed=True)
        self.overlap = Parameter.default(self.overlap, name=f'{self.name} overlap', fixed=True)

        super().__post_init__()

@dataclass
class SubstrateInterface(BaseGroupInterface):
    """Refl1D interface for Box2Err, specifically when used as a base group
    """

    _molgroup: Box2Err | None = None

    rho: Parameter = field(default_factory=lambda: Parameter(name='rho substrate', value=2.07))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='substrate roughness', value=2.07))

    substrate_surface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='substrate_surface', description='surface of substrate'))

    def __post_init__(self) -> None:
        self._molgroup = Box2Err(name=self.name)

        self.substrate_surface.set_function(functools.partial(lambda box: box.z + 0.5 * box.length, self._molgroup))

        super().__post_init__()
        
    def update(self, bulknsld: float):

        self._molgroup.fnSet(volume=self.normarea.value * self.overlap.value * 2.0,
                             length=2.0 * self.overlap.value,
                             position=0.0,
                             sigma=self.sigma.value,
                             nf=1.0,
                             nSL=self.normarea.value * self.overlap.value * 2.0 * self.rho.value)

@dataclass
class ssBLMInterface(BaseGroupInterface):
    """Refl1D interactor for ssBLM class
    """

    _molgroup: ssBLM | None = None

    lipids: List[Lipid] = field(default_factory=list)
    inner_lipid_nf: List[Parameter] = field(default_factory=list)
    outer_lipid_nf: List[Parameter] = field(default_factory=list)
    rho_substrate: Parameter = field(default_factory=lambda: Parameter(name='rho substrate', value=2.07))
    rho_siox: Parameter = field(default_factory=lambda: Parameter(name='rho siox', value=3.3))
    l_siox: Parameter = field(default_factory=lambda: Parameter(name='siox thickness', value=10))
    vf_bilayer: Parameter = field(default_factory=lambda: Parameter(name='volume fraction bilayer', value=0.9))
    l_lipid1: Parameter = field(default_factory=lambda: Parameter(name='inner acyl chain thickness', value=10.0))
    l_lipid2: Parameter = field(default_factory=lambda: Parameter(name='outer acyl chain thickness', value=10.0))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='bilayer roughness', value=5))
    global_rough: Parameter = field(default_factory=lambda: Parameter(name ='substrate roughness', value=5))
    l_submembrane: Parameter = field(default_factory=lambda: Parameter(name='submembrane thickness', value=10))

    substrate_surface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='substrate_surface', description='surface of substrate'))
    siox_surface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='siox_surface', description='surface of siox layer'))
    bilayer_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='bilayer_center', description='center of bilayer'))
    inner_headgroup_bottom: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_headgroup_bottom', description='bottom of inner headgroups'))
    inner_headgroup_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_headgroup_center', description='center of inner headgroups'))
    inner_hydrophobic_interface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_hydrophobic_interface', description='interface between inner headgroups and acyl chains'))
    outer_hydrophobic_interface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_hydrophobic_interface', description='interface between outer headgroups and acyl chains'))
    outer_headgroup_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_headgroup_center', description='center of outer headgroups'))
    outer_headgroup_top: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_headgroup_top', description='top of outer headgroups'))

    def __post_init__(self):
        self._molgroup = ssBLM(inner_lipids=self.lipids,
                               outer_lipids=self.lipids,
                             inner_lipid_nf=[p.value if hasattr(p, 'value') else p for p in self.inner_lipid_nf],
                             outer_lipid_nf=[p.value if hasattr(p, 'value') else p for p in self.outer_lipid_nf],
                             name=self.name)

        n_lipids = len(self.lipids)
        self._group_names = {'substrate': [f'{self.name}.substrate'],
                'silicon dioxide': [f'{self.name}.siox'],
                f'{self.name} inner headgroups': [f'{self.name}.headgroup1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} inner acyl chains': [f'{self.name}.methylene1_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer acyl chains': [f'{self.name}.methylene2_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl2_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer headgroups': [f'{self.name}.headgroup2_{i}' for i in range(1, n_lipids + 1)],
                }

        # connect reference points
        self.substrate_surface.set_function(functools.partial(lambda blm: blm.substrate.z + 0.5 * blm.substrate.length, self._molgroup))
        self.siox_surface.set_function(functools.partial(lambda blm: blm.siox.z + 0.5 * blm.siox.length, self._molgroup))
        self.bilayer_center.set_function(self._molgroup.fnGetCenter)
        self.inner_headgroup_bottom.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc - blm.av_hg1_l, self._molgroup))
        self.inner_headgroup_center.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc - 0.5 * blm.av_hg1_l, self._molgroup))
        self.inner_hydrophobic_interface.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc, self._molgroup))
        self.outer_hydrophobic_interface.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc, self._molgroup))
        self.outer_headgroup_center.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc + 0.5 * blm.av_hg2_l, self._molgroup))
        self.outer_headgroup_top.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc + blm.av_hg2_l, self._molgroup))

        super().__post_init__()

    def update(self, bulknsld: float):

        self._molgroup.substrate.length = 2.0 * self.overlap.value

        self._molgroup.fnSet(sigma=self.sigma.value,
            bulknsld=bulknsld * 1e-6,
            global_rough=self.global_rough.value,
            rho_substrate=self.rho_substrate.value * 1e-6,
            rho_siox=self.rho_siox.value * 1e-6,
            l_lipid1=self.l_lipid1.value,
            l_lipid2=self.l_lipid2.value,
            l_siox=self.l_siox.value,
            l_submembrane=self.l_submembrane.value,
            vf_bilayer=self.vf_bilayer.value,
            nf_inner_lipids=[p.value for p in self.inner_lipid_nf],
            nf_outer_lipids=[p.value for p in self.outer_lipid_nf],
            radius_defect=1e8)
        
        self.normarea.value = self._molgroup.normarea

@dataclass
class tBLMInterface(BaseGroupInterface):
    """Refl1D interactor for ssBLM class
    """

    _molgroup: tBLM | None = None

    tether: Tether = field(default_factory=Tether)
    filler: Component = field(default_factory=lambda: bme)
    lipids: List[Lipid] = field(default_factory=list)
    inner_lipid_nf: List[Parameter] = field(default_factory=list)
    outer_lipid_nf: List[Parameter] = field(default_factory=list)
    rho_substrate: Parameter = field(default_factory=lambda: Parameter(name='rho substrate', value=2.07))
    vf_bilayer: Parameter = field(default_factory=lambda: Parameter(name='volume fraction bilayer', value=0.9))
    l_lipid1: Parameter = field(default_factory=lambda: Parameter(name='inner acyl chain thickness', value=10.0))
    l_lipid2: Parameter = field(default_factory=lambda: Parameter(name='outer acyl chain thickness', value=10.0))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='bilayer roughness', value=5))
    substrate_rough: Parameter = field(default_factory=lambda: Parameter(name ='substrate roughness', value=5))
    l_tether: Parameter = field(default_factory=lambda: Parameter(name='tether length', value=10))
    nf_tether: Parameter = field(default_factory=Parameter(name='number fraction tether', value=0.45)) # number fraction of tether molecules in inner leaflet
    mult_tether: Parameter = field(default_factory=Parameter(name='bME to tether ratio', value=3)) #ratio of bME to tether molecules at surface

    substrate_surface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='substrate_surface', description='surface of substrate'))
    filler_surface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='filler_surface', description='surface of filler molecule'))
    bilayer_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='bilayer_center', description='center of bilayer'))
    inner_headgroup_bottom: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_headgroup_bottom', description='bottom of inner headgroups'))
    inner_headgroup_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_headgroup_center', description='center of inner headgroups'))
    inner_hydrophobic_interface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='inner_hydrophobic_interface', description='interface between inner headgroups and acyl chains'))
    outer_hydrophobic_interface: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_hydrophobic_interface', description='interface between outer headgroups and acyl chains'))
    outer_headgroup_center: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_headgroup_center', description='center of outer headgroups'))
    outer_headgroup_top: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='outer_headgroup_top', description='top of outer headgroups'))

    def __post_init__(self):
        self._molgroup = tBLM(tether=self.tether,
                              filler=self.filler,
                              inner_lipids=self.lipids,
                              outer_lipids=self.lipids,
                              inner_lipid_nf=[p.value if hasattr(p, 'value') else p for p in self.inner_lipid_nf],
                              outer_lipid_nf=[p.value if hasattr(p, 'value') else p for p in self.outer_lipid_nf],
                              name=self.name)

        n_lipids = len(self.lipids)
        self._group_names = {'substrate': [f'{self.name}.substrate'],
                f'{self.name} bME': [f'{self.name}.bME'],
                f'{self.name} tether': [f'{self.name}.tether_bme', f'{self.name}.tether_free', f'{self.name}.tether_hg'],
                f'{self.name} tether acyl chains': [f'{self.name}.tether_methylene', f'{self.name}.tether_methyl'],
                f'{self.name} inner headgroups': [f'{self.name}.headgroup1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} inner acyl chains': [f'{self.name}.methylene1_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer acyl chains': [f'{self.name}.methylene2_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl2_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer headgroups': [f'{self.name}.headgroup2_{i}' for i in range(1, n_lipids + 1)],
                }

        self.substrate_surface.set_function(functools.partial(lambda blm: blm.substrate.z + 0.5 * blm.substrate.length, self._molgroup))
        self.filler_surface.set_function(functools.partial(lambda blm: blm.bme.z + 0.5 * blm.bme.length, self._molgroup))
        self.bilayer_center.set_function(self._molgroup.fnGetCenter)
        self.inner_headgroup_bottom.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc - blm.av_hg1_l, self._molgroup))
        self.inner_headgroup_center.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc - 0.5 * blm.av_hg1_l, self._molgroup))
        self.inner_hydrophobic_interface.set_function(functools.partial(lambda blm: blm.z_ihc - 0.5 * blm.l_ihc, self._molgroup))
        self.outer_hydrophobic_interface.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc, self._molgroup))
        self.outer_headgroup_center.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc + 0.5 * blm.av_hg2_l, self._molgroup))
        self.outer_headgroup_top.set_function(functools.partial(lambda blm: blm.z_ohc + 0.5 * blm.l_ohc + blm.av_hg2_l, self._molgroup))    

        super().__post_init__()

    def update(self, bulknsld: float):

        self._molgroup.substrate.length = 2.0 * self.overlap.value

        self._molgroup.fnSet(sigma=self.sigma.value,
            bulknsld=bulknsld * 1e-6,
            global_rough=self.substrate_rough.value,
            rho_substrate=self.rho_substrate.value * 1e-6,
            l_lipid1=self.l_lipid1.value,
            l_lipid2=self.l_lipid2.value,
            l_tether=self.l_tether.value,
            vf_bilayer=self.vf_bilayer.value,
            nf_tether=self.nf_tether.value,
            mult_tether=self.mult_tether.value,
            nf_inner_lipids=[p.value for p in self.inner_lipid_nf],
            nf_outer_lipids=[p.value for p in self.outer_lipid_nf],
            radius_defect=1e8)
        
        self.normarea.value = self._molgroup.normarea

# ============= Box-type objects ===============

class ComponentBoxInterface(MolgroupsInterface):
    pass

class ProteinBoxInterface(MolgroupsInterface):
    pass

class TetheredBoxInterface(MolgroupsInterface):
    pass

class TetheredBoxDoubleInterface(MolgroupsInterface):
    pass

# ============= Spline objects ================
@dataclass
class BoxHermiteInterface(MolgroupsInterface):
    
    _molgroup: BoxHermite | None = None
    _bulknsld: float = 0.0

    dSpacing: float = 15.0
    startz: Parameter = field(default_factory=lambda: Parameter(name='start position', value=20))
    Dp: List[Parameter] = field(default_factory=[])
    Vf: List[Parameter] = field(default_factory=[])
    rhoH: Parameter = field(default_factory=lambda: Parameter(name='rhoH', value=0.0))
    rhoD: Parameter = field(default_factory=lambda: Parameter(name='rhoD', value=0.0))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='roughness', value=5))

    center_of_volume: ReferencePoint = field(default_factory=lambda: ReferencePoint(name='center_of_volume', description='center of volume'))
    rho: ReferencePoint = field(default_factory=lambda: ReferencePoint(name=f'nSLD', description='H/D aware nSLD of Hermite'))

    def __post_init__(self):
        self._molgroup = BoxHermite(name=self.name, n_box=21)
        self._group_names = {f'{self.name}': [f'{self.name}']}

        self.center_of_volume.set_function(self._center_of_volume)
        self.rho.set_function(functools.partial(lambda self: sld_from_bulk(self.rhoH.value, self.rhoD.value, self._bulknsld), self))

        super().__post_init__()

    def update(self, bulknsld: float) -> None:

        self._bulknsld = bulknsld

        self._molgroup.fnSetRelative(dSpacing=self.dSpacing,
                                     dStart=self.startz.value,
                                     dDp=[d.value for d in self.Dp],
                                     dVf=[d.value for d in self.Vf],
                                     dnSLD=self.rho.value * 1e-6,
                                     dnf=self.nf.value,
                                     sigma=self.sigma.value)
        
# ============= Euler objects =================

class ContinuousEulerInterface(MolgroupsInterface):
    pass

# ============= Complex objects ===============

@dataclass
class BLMProteinComplexInterface(BaseGroupInterface):
    
    _molgroup: BLMProteinComplex | None = None

    base_blm: ssBLMInterface | tBLMInterface | None = None
    blms: List[BLMInterface] = field(default_factory=list)
    proteins: List[ProteinBoxInterface | BoxHermiteInterface] = field(default_factory=list)

    def __post_init__(self) -> None:

        self._molgroup = BLMProteinComplex(blms=[blm._molgroup for blm in self.all_blms],
                                           proteins=[prot._molgroup for prot in self.proteins])
        
        # compile group names
        _group_names = {}
        for prepend, gplist in zip(['blms', 'proteins'], [self.all_blms, self.proteins]):
            prepend = f'{self.name}.{prepend}'
            for gp in gplist:
                for k, gpnames in gp._group_names.items():
                    gpnames = [f'{prepend}.{gpname}' for gpname in gpnames]
                    _group_names.update({k: gpnames})
        self._group_names = _group_names

        print(self._group_names)
        #for gp in self.all_blms + self.proteins:
        #    self._group_names.update(gp._group_names)

        super().__post_init__()

        # tie base group overlap to this overlap, after conversion to a parameter
        self.base_blm.overlap = self.overlap

    def _get_parameters(self) -> Dict[str, Parameter]:

        pars = {}
        for gp in self.all_blms + self.proteins:
            pars.update(gp._get_parameters())
        return pars

    @property
    def all_blms(self) -> List[BLMInterface | ssBLMInterface | tBLMInterface]:
        return [self.base_blm] + self.blms if self.base_blm is not None else self.blms

    def update(self, bulknsld: float) -> None:

        for gp in self.all_blms + self.proteins:
            gp.update(bulknsld)

        self._molgroup.fnAdjustBLMs()
        self.normarea.value = self._molgroup.normarea

    def store_profile(self, z: np.ndarray) -> Dict:
        # special profile storage that takes into account excess density.
        # TODO: this is somewhat hackish. It might make more sense to have a subclass of MolgroupsLayer that
        # incorporates this logic
        super().store_profile(z)

        normarea = self.normarea.value
        
        prot_area = np.zeros_like(z)
        prot_nsl = np.zeros_like(z)
        for gp in self.proteins:
            _, area, nsl = gp.render(z)
            prot_area += area
            prot_nsl += nsl
        
        blm_area = np.zeros_like(z)
        blm_nsl = np.zeros_like(z)
        for gp in self.all_blms:
            _, area, nsl = gp.render(z)
            blm_area += area
            blm_nsl += nsl
        
        frac_replacement = np.ones_like(area)
        if len(self.proteins):
            over_filled = (blm_area + prot_area) > normarea
            frac_replacement[over_filled] = (blm_area / (normarea - prot_area))[over_filled]

        for blm in self.all_blms:
            for gplist in blm._group_names.values():
                for gp in gplist:
                    self._stored_profile[f'{self.name}.blms.{gp}']['area'] /= frac_replacement
                    self._stored_profile[f'{self.name}.blms.{gp}']['sl'] /= frac_replacement

        self._stored_profile['area'] = blm_area / frac_replacement + prot_area
        self._stored_profile['sl'] = blm_nsl / frac_replacement + prot_nsl
        self._stored_profile['normarea'] = normarea
