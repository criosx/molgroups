"""Module for interfacing with Refl1D
    https://github.com/reflectometry/refl1d
"""

from typing import ClassVar, List, Tuple
from dataclasses import dataclass, InitVar, field, fields, asdict
import uuid

import numpy as np

from refl1d.names import Slab, Stack, Parameter, Experiment, FitProblem, SLD
from refl1d.flayer import FunctionalProfile
from refl1d.model import Layer
from refl1d.probe import QProbe

from molgroups.mol import nSLDObj, ssBLM, tBLM
from molgroups.components import Component, Lipid, Tether, bme

def write_groups(groups: List[nSLDObj], labels: List[str], dimension: int, stepsize: float):
    """Return dictionaries with combined output of fnWriteGroup2Dict and fnWriteResults2Dict
    
        Inputs:
        groups: list of Molgroups objects to process
        labels: list (same length as groups) of labels"""
    
    moldict = {}
    resdict = {}
    for lbl, gp in zip(labels, groups):
        moldict = {**moldict, **gp.fnWriteGroup2Dict({}, lbl, np.arange(dimension) * stepsize)}
        resdict = {**resdict, **gp.fnWriteResults2Dict({}, lbl)}
        
    return moldict, resdict

@dataclass
class BaseMolgroupsInterface:
    """Base class for interacting with molgroups objects
    """

    id: str | None = None
    name: str | None = None
    _molgroup: nSLDObj | None = None
    _stored_profile: dict | None = None
    _group_names: dict[str, List[str]] = field(default_factory=dict)

    def __post_init__(self) -> None:

        if self.id is None:
            self.id = str(uuid.uuid4())

        for p in self._get_parameters().values():
            if not p.name.startswith(self.name):
                p.name = self.name + ' ' + p.name

    def _get_parameters(self) -> dict[str, Parameter]:
        """Gets a list of the parameters associated with the interactor

        Returns:
            List[Parameter]: Parameter list
        """
        parlist = [f.name for f in fields(self) if f.type == Parameter]
        return {p: getattr(self, p) for p in parlist}

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

@dataclass
class ssBLMInterface(BaseMolgroupsInterface):
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

@dataclass
class tBLMInterface(BaseMolgroupsInterface):
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
        

# =========================

@dataclass(init=False)
class MolgroupsLayer(Layer):

    base_group: BaseMolgroupsInterface
    add_groups: List[BaseMolgroupsInterface] = field(default_factory=list)
    overlay_groups: List[BaseMolgroupsInterface] = field(default_factory=list)
    contrast: SLD | None=None
    overlap: float | Parameter = 20.0
    thickness: float | Parameter = 0.0
    name: str | None = None

    def __init__(self,
                 base_group: BaseMolgroupsInterface,
                 add_groups: List[BaseMolgroupsInterface] = [],
                 overlay_groups: List[BaseMolgroupsInterface] = [],
                 contrast: SLD | None=None,
                 overlap: float | Parameter = 20.0,
                 thickness: float | Parameter = 0.0,
                 name=None) -> None:

        self.base_group = base_group
        self.add_groups = add_groups
        self.overlay_groups = overlay_groups
        self.overlap = overlap

        if name is None:
            name = self.base_group.name
        self.overlap = Parameter.default(overlap, name=name+ " overlap")
        self.thickness = Parameter.default(thickness, name=name+" thickness")
        self.interface = Parameter.default(0.0, name=name+" interface")

        self.name = name
        self.magnetism = None
        self.contrast = contrast

    def update(self):

        # 1. update base_group
        self.base_group.update(bulknsld=self.contrast.rho.value)
        normarea = self.base_group._molgroup.normarea \
                    if hasattr(self.base_group._molgroup, 'normarea') \
                        else 1.0

        # 2. apply normarea to all remaining objects
        for group in self.add_groups + self.overlay_groups:
            if hasattr(group._molgroup, 'normarea'):
                group._molgroup.fnSetNormarea(normarea)
        
        # 3. update all remaining objects
        for group in self.add_groups + self.overlay_groups:
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

    def _filled_profile(self, z):
        """Given area and nSL profiles, fill in the remaining volume with bulk material"""
        
        self.update()
        normarea, area, nsl = self.profile(z)

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

    def _custom_plot(self):

        import plotly.graph_objs as go
        from refl1d.webview.server.colors import COLORS

        def hex_to_rgb(hex_string):
            r_hex = hex_string[1:3]
            g_hex = hex_string[3:5]
            b_hex = hex_string[5:7]
            return int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)

        # compile moldat
        moldat = {}
        group_names = {}
        for group in [self.base_group] + self.add_groups + self.overlay_groups:
            for k, v in group._group_names.items():
                # propagate frac_replacement to subgroups
                for gp in v:
                    group._stored_profile[gp]['frac_replacement'] = group._stored_profile['frac_replacement']

                # merge group names
                if k not in group_names.keys():
                    group_names[k] = []
                group_names[k] += v

            
            moldat.update(group._stored_profile)

        # define normarea
        normarea = self.base_group._stored_profile['normarea']

        fig = go.Figure()
        traces = []
        MOD_COLORS = COLORS[1:]
        color_idx = 1
        sumarea = 0
        for lbl, item in group_names.items():
            area = 0
            for gp in item:
                if gp in moldat.keys():
                    zaxis = moldat[gp]['zaxis']
                    newarea = moldat[gp]['area'] / moldat[gp]['frac_replacement']
                    area += np.maximum(0, newarea)
                else:
                    print(f'Warning: {gp} not found')

            color = MOD_COLORS[color_idx % len(MOD_COLORS)]
            plotly_color = ','.join(map(str, hex_to_rgb(color)))
            traces.append(go.Scatter(x=zaxis,
                                    y=area / normarea,
                                    mode='lines',
                                    name=lbl,
                                    line=dict(color=color)))
            traces.append(go.Scatter(x=zaxis,
                                    y=area / normarea,
                                    mode='lines',
                                    line=dict(width=0),
                                    fill='tozeroy',
                                    fillcolor=f'rgba({plotly_color},0.3)',
                                    showlegend=False
                                    ))
            color_idx += 1
            sumarea += area

        color = COLORS[0]
        plotly_color = ','.join(map(str, hex_to_rgb(color)))
        
        traces.append(go.Scatter(x=zaxis,
                                    y=sumarea / normarea,
                                    mode='lines',
                                    name='buffer',
                                    line=dict(color=color)))
        traces.append(go.Scatter(x=zaxis,
                                    y=sumarea / normarea,
                                    mode='lines',
                                    line=dict(width=0),
                                    fill='tonexty',
                                    fillcolor=f'rgba({plotly_color},0.3)',
                                    showlegend=False
                                    ))    
        traces.append(go.Scatter(x=zaxis,
                                    y=[1.0] * len(zaxis),
                                    mode='lines',
                                    line=dict(color=color, width=0),
                                    showlegend=False))

        
        fig.add_traces(traces[::-1])

        fig.update_layout(
            title='Component Volume Occupancy',
            template = 'plotly_white',
            xaxis_title=dict(text='z (Ang)'),
            yaxis_title=dict(text='volume occupancy')
        )

        return fig

# =============

@dataclass(init=False, eq=False)
class MolgroupsStack(Stack):

    substrate: Stack | Slab
    molgroups_layer: MolgroupsLayer

    def __init__(self,
                 substrate: Stack | Slab,
                 molgroups_layer: MolgroupsLayer,
                 name="MolgroupsStack",
                 **kw):
        self.substrate, self.molgroups_layer = substrate, molgroups_layer
        layer_contrast = Slab(material=molgroups_layer.contrast, thickness=0.0000, interface=0.0000)
        if isinstance(substrate, Stack):
            substrate.layers[-1].thickness = substrate.layers[-1].thickness - molgroups_layer.overlap
        elif isinstance(substrate, Slab):
            substrate.thickness = substrate.thickness - molgroups_layer.overlap
        layers = substrate | molgroups_layer | layer_contrast
        super().__init__(base=None, 
                         layers=layers,
                         name=name,
                         interface=None,
                         thickness=None)

# =============

@dataclass(init=False)
class MolgroupsExperiment(Experiment):

    sample: MolgroupsStack | None = None,

    def __init__(self,
                 sample: MolgroupsStack | None = None,
                 probe=None,
                 name=None,
                 roughness_limit=0,
                 dz=None,
                 dA=None,
                 step_interfaces=None,
                 smoothness=None,
                 interpolation=0,
                 constraints=None,
                 version: str | None = None,
                 auto_tag=False):
        super().__init__(sample, probe, name, roughness_limit, dz, dA, step_interfaces, smoothness, interpolation, constraints, version, auto_tag)
        self._custom_plot = self.sample.molgroups_layer._custom_plot

def make_samples(layer_template: MolgroupsLayer, substrate: Stack, contrasts: List[SLD]) -> List[MolgroupsExperiment]:
    """Create samples from combining a substrate stack with a molgroups layer
    
        Args:
            layer_template: molgroups layer template
            substrate: Refl1D Stack or Layer object representing the substrate
            contrasts: list of buffer materials, e.g. [d2o, h2o]. One sample will be created for each contrast
    """
    samples = []

    for contrast in contrasts:
        mollayer = MolgroupsLayer(base_group=layer_template.base_group,
                                  add_groups=layer_template.add_groups,
                                  overlay_groups=layer_template.overlay_groups,
                                  contrast=contrast,
                                  overlap=layer_template.overlap,
                                  thickness=layer_template.thickness,
                                  name=contrast.name + ' ' + layer_template.name)
        samples.append(MolgroupsStack(substrate=substrate,
                                      molgroups_layer=mollayer))

    return samples

if __name__ == '__main__':
    from bumps.serialize import serialize, deserialize
    from molgroups.lipids import DOPC, DOPE
    import matplotlib.pyplot as plt

    vf_bilayer = Parameter(name='volume fraction bilayer', value=0.9).range(0.0, 1.0)
    vf_bilayer_peptide = Parameter(name='volume fraction bilayer peptide', value=0.9).range(0.0, 1.0)
    l_lipid1 = Parameter(name='inner acyl chain thickness', value=10.0).range(8, 30)
    l_lipid2 = Parameter(name='outer acyl chain thickness', value=10.0).range(8, 18)
    dl_lipid_peptide = Parameter(name='acyl chain thickness difference peptide', value=0.0).range(-3, 3)
    dl_lipid_low = Parameter(name='acyl chain thickness difference low salt', value=0.0) #.range(-3, 3)
    #sigma = Parameter(name='bilayer roughness', value=5).range(0.5, 9)
    global_rough = Parameter(name ='substrate roughness', value=5).range(0.5, 9)
    sigma = global_rough

    d_oxide = Parameter(name='silicon oxide layer thickness', value=10).range(5, 30)
    d_Cr =  Parameter(name='chromium layer thickness', value=40).range(10, 150)
    d_gold =  Parameter(name='gold layer thickness', value=140).range(100, 200) #thickness of gold
    rough_cr_au =  Parameter(name='gold chromium roughness', value=10).range(2, 24.0) # roughness of Cr/Au interface
    l_submembrane = Parameter(name='submembrane thickness', value=10).range(0, 50)
    l_submembrane_low = Parameter(name='submembrane thickness low salt', value=10).range(0, 50)
    l_submembrane_peptide = Parameter(name='submembrane thickness peptide', value=10).range(0, 50)

    ## Substrate materials
    silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)

    siox = SLD(name='siox', rho=4.1000, irho=0.0000)
    siox.rho.range(2.9000, 5.1000)

    cr = SLD(name='chromium', rho=2.7, irho=0.0)
    cr.rho.range(2.2000, 4.0000)

    gold = SLD(name='gold', rho=4.4, irho=0.0) #iro is the absorption of neutrons, should be 0
    gold.rho.range(4.2000, 4.8000)

    gasket = SLD(name='gasket', rho=0, irho=0.0)
    gasket.rho.range(-0.5, 6.4)

    d2o = SLD(name='d2o', rho=6.3000, irho=0.0000)
    d2o_low = SLD(name='d2o low salt', rho=6.3000, irho=0.0000)
    d2o_peptide = SLD(name='d2o peptide', rho=6.3000, irho=0.0000)
    h2o = SLD(name='h2o', rho=-0.56, irho=0.0000)
    h2o_peptide = SLD(name='h2o peptide', rho=-0.56, irho=0.0000)

    # nSLD parameters
    d2o.rho.range(6.0000, 6.5000)
    d2o_low.rho.range(6.0000, 6.5000)
    d2o_peptide.rho.range(5.6, 6.5)
    h2o.rho.range(-0.56, 2.6)
    h2o_peptide.rho.range(-0.56, 2.6)

    ## Then bulk layers are created, each with its own 'material'.  If you want to force
    ## two layers to always match SLD you can use the same material in multiple layers.
    ## The roughnesses of each layer are set to zero to begin with:

    layer_silicon = Slab(material=silicon, thickness=0.0000, interface=global_rough)
    layer_siox = Slab(material=siox, thickness=d_oxide, interface=global_rough)
    layer_cr = Slab(material=cr, thickness=d_Cr, interface=rough_cr_au)
    layer_gold = Slab(material=gold, thickness=d_gold, interface=0.0000)
    layer_gasket = Slab(material=gasket, thickness=0, interface=global_rough)

    substrate = layer_silicon | layer_siox | layer_cr | layer_gold

    blm = ssBLMInterface(name='bilayer',
                          lipids=[DOPC, DOPE],
                          lipid_nf=[0.6, 0.4],
                          l_submembrane=l_submembrane)
    
    blm.rho_substrate = gold.rho
        
    mlayer = MolgroupsLayer(base_group=blm,
                            contrast=d2o,
                            thickness=100)
    
    #print(mlayer.parameters())
    
    mstack = MolgroupsStack(substrate, mlayer)
    
    Q = np.linspace(0.01, 0.1, 101)
    dQ = 0.02 * Q
    mexp = MolgroupsExperiment(sample=mstack, probe=QProbe(Q, dQ), dz=1.0)

    problem=FitProblem(mexp)

    setattr(problem, 'custom_plot', mexp._custom_plot)
    #print(problem.labels())

    #plt.figure()
    #problem.plot()
    #plt.show()
