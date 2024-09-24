"""Plots for the Refl1D webview interface
    Used with MolgroupsExperiment.register_webview_plot
"""
import numpy as np
import plotly.graph_objs as go

from bumps.dream.state import MCMCDraw
from bumps.webview.server.custom_plot import CustomWebviewPlot
from refl1d.names import FitProblem, Experiment
from refl1d.webview.server.colors import COLORS

def hex_to_rgb(hex_string):
    r_hex = hex_string[1:3]
    g_hex = hex_string[3:5]
    b_hex = hex_string[5:7]
    return int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)

from .layers import MolgroupsLayer

def cvo_plot(layer: MolgroupsLayer, model: Experiment | None = None, problem: FitProblem | None = None, state: MCMCDraw | None = None):
    # component volume occupancy plot

    # compile moldat
    moldat = {}
    group_names = {}
    for group in [layer.base_group] + layer.add_groups + layer.overlay_groups:
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
    normarea = layer.base_group._stored_profile['normarea']

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
        yaxis_title=dict(text='volume occupancy'),
        legend=dict(yanchor='top',
                    xanchor='right',
                    x=0.99,
                    y=0.99),
        yaxis_range=[0, 1]
    )

    return CustomWebviewPlot(fig_type='plotly',
                                plotdata=fig)