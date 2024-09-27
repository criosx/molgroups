"""Support for directly interfacing the base molgroups objects via
    FunctionalProfile.
   Notes:
    1. A FitProblem defined this way is not serializable and does
        not automatically implement the webview custom plots
    2. Currently the register_cvo_plot does not work perfectly.
        The molgroups objects are attached to the model but are then
        copied before plotting. So the plot objects don't point back
        to the original object and the plot ends up being one plot
        refresh behind. This can be seen when switching between different
        CVO plots. The underlying molgroups objects are not updated
        in the copied models.
"""
from typing import List, Dict, Callable, Tuple

import numpy as np

from ..mol import nSLDObj

from bumps.dream.state import MCMCDraw
from bumps.webview.server.custom_plot import CustomWebviewPlot
from refl1d.names import Slab, Stack, SLD, Experiment, FitProblem
from refl1d.flayer import FunctionalProfile
from refl1d.webview.server.colors import COLORS

import plotly.graph_objs as go

def hex_to_rgb(hex_string):
    r_hex = hex_string[1:3]
    g_hex = hex_string[3:5]
    b_hex = hex_string[5:7]
    return int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)

def write_groups(z: list | np.ndarray,
                 groups: List[nSLDObj],
                 labels: List[str]) -> Tuple[Dict, Dict]:
    """Return dictionaries with combined output of fnWriteGroup2Dict and fnWriteResults2Dict
    
    Args:
        z (List | np.ndarray): list or np.array of z values
        groups (List[nSLDObj]): list of Molgroups objects to process
        labels (List[str]): list (same length as groups) of labels

    Returns:
        Dict: output of nSLDObj.fnWriteGroup2Dict, keys are object
            names
        Dict: output of nSLDobj.fnWriteResults2Dict, keys are object
            names
    """
    
    moldict = {}
    resdict = {}
    for lbl, gp in zip(labels, groups):
        moldict = {**moldict, **gp.fnWriteGroup2Dict({}, lbl, z)}
        resdict = {**resdict, **gp.fnWriteResults2Dict({}, lbl)}
        
    return moldict, resdict

def apply_bulknsld(z: np.ndarray,
                   bulknsld: float,
                   normarea: float | np.ndarray,
                   area: np.ndarray,
                   nsl: np.ndarray) -> np.ndarray:
    
    """Given area and nSL profiles, fill in the remaining volume with bulk material.

    Args:
        z (np.ndarray): spatial domain over which to render molecular
            groups
        bulknsld (float): scattering length density of medium
        normarea (float | np.ndarray): area normalization factor. If 
            array, should be a constant value
        area (np.ndarray): area taken up by molecular groups. The total
            volume fraction is area / normarea
        nsl (np.ndarray): total scattering length of molecular groups,
            defined as nSLD * area.

    Returns:
        np.ndarray: nSLD profile
    """

    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * np.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # Return nSLD profile in Refl1D units
    return nsld

def make_samples(func: Callable,
                 npoints: int,
                 substrate: Stack | Slab,
                 contrasts: List[SLD],
                 **kwargs):
    """Create samples from combining a substrate stack with a molgroups
        layer capped with a bulk contrast layer
    
    Args:
        func (Callable): function used to define FunctionalProfile
            object. Must have form func(z, bulknsld, *args)
        npoints (int): number of points in FunctionalProfile.
            Equal to desired layer thickness / (Experiment dz)
        substrate (Stack | Slab): Refl1D Stack or Slab object
            representing the layer below the lowest molecular group
        contrasts (List[SLD]): list of buffer materials, e.g.
            [d2o, h2o]. One sample will be created for each contrast
        **kwargs: keyword arguments. Must have one keyword argument for
            each arg in func(..., *args), but not one for bulknsld
        
    Returns:
        List[Stack]: List of Refl1D samples, one for each contrast.
            Stack has structure 'substrate | mollayer | contrast' so
            the last group is a semiinfinite layer of the contrast
            medium
    """
            
    samples = []

    for contrast in contrasts:
        mollayer = FunctionalProfile(npoints, 0, profile=func, bulknsld=contrast.rho, **kwargs)
        layer_contrast = Slab(material=contrast, thickness=0.0000, interface=5.0000)
        samples.append(substrate | mollayer | layer_contrast)

    return samples

class CVOPlot:

    def __init__(self,
                 z: np.ndarray,
                 groups: List[nSLDObj],
                 labels: List[str],
                 group_names: Dict[str, List[str]],
                 normarea_group: str | None = None) -> None:
    
        self.z = z
        self.groups = groups
        self.labels = labels
        self.group_names = group_names
        self.normarea_group = normarea_group

    def __call__(self, model, problem, state):
        return cvo_functionalprofileplot(self.z, self.groups, self.labels,self.group_names, model, problem, state, self.normarea_group)

def register_cvo_plot(model: Experiment,
                  z: np.ndarray,
                  groups: List[nSLDObj],
                  labels: List[str],
                  group_names: Dict[str, List[str]],
                  normarea_group: str | None = None) -> None:
    """Attaches a CVO plot to a Refl1D Experiment

    Args:
        model (Experiment): Refl1D Experiment to which to attach the
            plot.
        z (np.ndarray): z-axis of molgroup layer
        groups (List[nSLDObj]): list of Molgroups objects to process
        labels (List[str]): list (same length as groups) of labels
        group_names (Dict[str, List[str]]): group names for plotting
            Each plotted group will be labeled with the key, and the sum
            of the areas of the groups in the list will be plotted
        normarea_group (str | None): Optional, default None.
            If None, defaults to 1, otherwise
            use defined group name, e.g. bilayer.normarea,
            for normalizing area.
    """

    # NOTE: must use lambda here, not functools.partial or CVOPlot
    model.register_webview_plot('Component Volume Occupancy',
                                lambda model, problem, state: cvo_functionalprofileplot(
                                        z,
                                        groups,
                                        labels,
                                        group_names,
                                        model, problem, state,
                                        normarea_group=normarea_group),
                                'parameter')

def cvo_functionalprofileplot(z: np.ndarray,
                              groups: List[nSLDObj],
                              labels: List[str],
                              group_names: Dict[str, List[str]],
                              model: Experiment,
                              problem: FitProblem,
                              state: MCMCDraw | None,
                              normarea_group: str | None = None,
                              ) -> CustomWebviewPlot:

    """Component volume occupancy profile plot for Refl1D webview GUI.
        
    Usage:
        from functools import partial
        model.register_webview_plot('<Plot Title>',
                                    partial(cvo_functionprofileplot,
                                            np.arange(DIMENSION) * STEPSIZE,
                                            [blm, other_group],
                                            ['bilayer', 'other_thing'],
                                            {'substrate': 'bilayer.substrate'},
                                            normarea_group='bilayer.normarea'),
                                    'parameter')

    Args:
        z (np.ndarray): z-axis of molgroup layer
        groups (List[nSLDObj]): list of Molgroups objects to process
        labels (List[str]): list (same length as groups) of labels
        group_names (Dict[str, List[str]]): group names for plotting
            Each plotted group will be labeled with the key, and the sum
            of the areas of the groups in the list will be plotted
        model (Experiment): Refl1D Experiment object. Provided by plotter.
        problem (FitProblem): Refl1D FitProblem object. Provided by plotter.
        state (bumps.dream.state.MCMCDraw): Uncertainty state. Provided
            by plotter.
        normarea_group (str | None): Optional, default None.
            If None, defaults to 1, otherwise
            use defined group name, e.g. bilayer.normarea,
            for normalizing area.

    Returns:
        CustomWebviewPlot: plotly plot object container interpretable
            by the Refl1D webview GUI
    """

    moldat, _ = write_groups(z, groups, labels)

    if normarea_group is None:
        normarea = 1.0
    else:
        normarea = moldat[normarea_group]['area']

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
                area += np.maximum(0, moldat[gp]['area'])
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

    
    fig.add_traces((traces)[::-1])

    fig.update_layout(
        title='Component Volume Occupancy',
        template = 'plotly_white',
        xaxis_title=dict(text='z (Ang)'),
        yaxis_title=dict(text='volume occupancy'),
                legend=dict(yanchor='top',
                    y=0.98,
                    xanchor='right',
                    x=0.99),
        yaxis_range=[0, 1]
    )

    return CustomWebviewPlot(fig_type='plotly',
                             plotdata=fig)