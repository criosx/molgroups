"""Refl1D shims for molgroups.mol objects.
    Uses Refl1D Parameter objects for parameters, thus allowing model serialization
"""

from typing import List
from dataclasses import dataclass, field, fields
import uuid

import numpy as np

from refl1d.names import Parameter

from molgroups.mol import nSLDObj, ssBLM, tBLM, Box2Err, BoxHermite
from molgroups.components import Component, Lipid, Tether, bme

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

# ============= BaseGroup objects ===============

@dataclass
class BaseGroupInterface(MolgroupsInterface):
    """Interface specifically for base groups, i.e. those that occupy the edges of the molgroups canvas
    """

    normarea: Parameter | float = 1.0
    overlap: Parameter | float = 20.0

    def __post_init__(self) -> None:

        self.normarea = Parameter.default(self.normarea, name=f'{self.name} normarea')
        self.overlap = Parameter.default(self.overlap, name=f'{self.name} overlap')

        super().__post_init__()

@dataclass
class SubstrateInterface(BaseGroupInterface):
    """Refl1D interface for Box2Err, specifically when used as a base group
    """

    _molgroup: Box2Err | None = None

    rho: Parameter = field(default_factory=lambda: Parameter(name='rho substrate', value=2.07))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='substrate roughness', value=2.07))

    def __post_init__(self) -> None:
        super().__post_init__()
        self._molgroup = Box2Err(name=self.name)
        
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
    lipid_nf: List[float] = field(default_factory=list)
    rho_substrate: Parameter = field(default_factory=lambda: Parameter(name='rho substrate', value=2.07))
    rho_siox: Parameter = field(default_factory=lambda: Parameter(name='rho siox', value=3.3))
    l_siox: Parameter = field(default_factory=lambda: Parameter(name='siox thickness', value=10))
    vf_bilayer: Parameter = field(default_factory=lambda: Parameter(name='volume fraction bilayer', value=0.9))
    l_lipid1: Parameter = field(default_factory=lambda: Parameter(name='inner acyl chain thickness', value=10.0))
    l_lipid2: Parameter = field(default_factory=lambda: Parameter(name='outer acyl chain thickness', value=10.0))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='bilayer roughness', value=5))
    global_rough: Parameter = field(default_factory=lambda: Parameter(name ='substrate roughness', value=5))
    l_submembrane: Parameter = field(default_factory=lambda: Parameter(name='submembrane thickness', value=10))

    def __post_init__(self):
        super().__post_init__()
        self._molgroup = ssBLM(lipids=self.lipids,
                              lipid_nf=self.lipid_nf)

        n_lipids = len(self.lipids)
        self._group_names = {'substrate': [f'{self.name}.substrate'],
                'silicon dioxide': [f'{self.name}.siox'],
                f'{self.name} inner headgroups': [f'{self.name}.headgroup1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} inner acyl chains': [f'{self.name}.methylene1_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl1_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer acyl chains': [f'{self.name}.methylene2_{i}' for i in range(1, n_lipids + 1)] + [f'{self.name}.methyl2_{i}' for i in range(1, n_lipids + 1)],
                f'{self.name} outer headgroups': [f'{self.name}.headgroup2_{i}' for i in range(1, n_lipids + 1)],
                }
        
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
    lipid_nf: List[float] = field(default_factory=list)
    rho_substrate: Parameter = field(default_factory=lambda: Parameter(name='rho substrate', value=2.07))
    vf_bilayer: Parameter = field(default_factory=lambda: Parameter(name='volume fraction bilayer', value=0.9))
    l_lipid1: Parameter = field(default_factory=lambda: Parameter(name='inner acyl chain thickness', value=10.0))
    l_lipid2: Parameter = field(default_factory=lambda: Parameter(name='outer acyl chain thickness', value=10.0))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='bilayer roughness', value=5))
    substrate_rough: Parameter = field(default_factory=lambda: Parameter(name ='substrate roughness', value=5))
    l_tether: Parameter = field(default_factory=lambda: Parameter(name='tether length', value=10))
    nf_tether: Parameter = field(default_factory=Parameter(name='number fraction tether', value=0.45)) # number fraction of tether molecules in inner leaflet
    mult_tether: Parameter = field(default_factory=Parameter(name='bME to tether ratio', value=3)) #ratio of bME to tether molecules at surface


    def __post_init__(self):
        super().__post_init__()
        self._molgroup = tBLM(tether=self.tether,
                              filler=self.filler,
                              lipids=self.lipids,
                              lipid_nf=self.lipid_nf)

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
            radius_defect=1e8)
        
        self.normarea.value = self._molgroup.normarea

class BaseBLMProteinComplexInterface(BaseGroupInterface):
    pass

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

    dSpacing: float = 15.0
    startz: Parameter = field(default_factory=lambda: Parameter(name='start position', value=20))
    Dp: List[Parameter] = field(default_factory=[])
    Vf: List[Parameter] = field(default_factory=[])
    rho: Parameter = field(default_factory=lambda: Parameter(name='rho', value=0.0))
    sigma: Parameter = field(default_factory=lambda: Parameter(name='roughness', value=5))

    def __post_init__(self):
        self._molgroup = BoxHermite(name=self.name, n_box=21)
        self._group_names = {f'{self.name}': [f'{self.name}']}

        super().__post_init__()

    def update(self, bulknsld: float) -> None:

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

class BaseBLMProteinComplexInterface(BaseGroupInterface):
    pass
