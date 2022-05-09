# -*- coding: utf-8 -*-
import plotly.graph_objects as go
import plotly.express.colors as pxc

my_template = go.layout.Template()

my_template.layout.autosize = False

my_template.layout.xaxis = {
    'title': {
        'font': {
            'color': 'black',
            'size': 22,
            'family': 'Avenir'
        },
        'standoff': 12
    },
    'tickfont': {
        'color': 'black',
        'size': 18
    },
    'automargin': True,
    'showgrid': False,
    'showline': True,
    'ticks': 'outside',
    'tickformat': '.1f'
}

my_template.layout.yaxis = {
    'title': {
        'font': {
            'color': 'black',
            'size': 22,
            'family': 'Avenir'
        },
        'standoff': 12
    },
    'tickfont': {
        'color': 'black',
        'size': 18,
        'family': 'Avenir'
    },
    'automargin': True,
    'showgrid': False,
    'showline': True,
    'ticks': 'outside',
    'tickformat': '.2f'
}

my_template.layout.legend = {
    'font': {
        'color': 'black',
        'size': 18,
        'family': 'Avenir'
    }
}

line_styles = {
    0: {'color': pxc.sequential.Sunsetdark[4], 'dash': 'solid'},
    0.25: {'color': pxc.qualitative.Pastel2[1], 'dash': 'dash'},
    0.5: {'color': pxc.qualitative.Pastel2[1], 'dash': 'dot'},
    0.75: {'color': pxc.qualitative.Pastel2[1], 'dash': 'dashdot'},
    1.0: {'color': pxc.qualitative.Set2[1], 'dash': 'solid'},
    1.25: {'color': pxc.qualitative.Pastel2[0], 'dash': 'dash'},
    1.50: {'color': pxc.qualitative.Pastel2[0], 'dash': 'dot'},
    1.75: {'color': pxc.qualitative.Pastel2[0], 'dash': 'dashdot'},
    2.0: {'color': pxc.qualitative.Dark2[0], 'dash': 'solid'},

    -0.25: {'color': pxc.qualitative.Pastel1[1], 'dash': 'dash'},
    -0.5: {'color': pxc.qualitative.Pastel1[1], 'dash': 'dot'},
    -0.75: {'color': pxc.qualitative.Pastel1[1], 'dash': 'dashdot'},
    -1.0: {'color': pxc.qualitative.Safe[0], 'dash': 'solid'},
    -1.25: {'color': pxc.qualitative.Dark2[-1], 'dash': 'dash'},
    -1.50: {'color': pxc.qualitative.Dark2[-1], 'dash': 'dot'},
    -1.75: {'color': pxc.qualitative.Dark2[-1], 'dash': 'dashdot'},
    -2.0: {'color': pxc.sequential.Sunsetdark[-1], 'dash': 'solid'},

    'orbital_zeeman': {'color': 'teal'},
    'spin_zeeman': {'color': 'brown'},
    'diamagnetic': {'color': 'gray'}
}


def visualize_mu_scan(df, df_Sz=None):
    """Create a figure for the 'onion' plots, i.e. (mu, R, S_z) after a mu PES scan."""


    figure = go.Figure()

    figure.update_layout(
        width=800, height=800,
        margin=dict(t=10, r=0, l=0, b=0)
    )


    # Add the (mu, R, S_z) scatter.
    R_values = df['R'].unique()
    for R in R_values:
        df_R = df[abs(df['R'] - R) < 1.0e-08]
        figure.add_trace(
            go.Scatter3d(
                x = df_R['mu'],
                y = df_R['R'],
                z = df_R['S_z'],
                mode = 'markers+lines',
                marker = {
                    'size': 2,
                    'color': df_R['S_z'],
                    'colorscale': 'Sunsetdark'
                },
                line = {
                    'color': 'black',
                    'width': 0.2
                },
                showlegend = False
            )
        )


    if df_Sz is not None:
        # Add the GHF Sz trace.
        figure.add_trace(
            go.Scatter3d(
                x = [0] * len(df_Sz['R']),
                y = df_Sz['R'],
                z = df_Sz['S_z'],
                mode = 'lines',
                line = {
                    'color': 'royalblue',
                    'width': 4
                },
                showlegend = False
            )
        )


    # Add axis titles.
    figure.update_scenes(
        xaxis_title = "μ",
        yaxis_title = "Internuclear distance (a.u.)",
        zaxis_title = "S<sub>z</sub> (a.u.)"
    )


    # Adjust the camera position.
    figure.update_layout(
        scene_camera = {
            'up': {'x': 0, 'y': 0, 'z': 1},
            'center': {'x': 0, 'y': 0, 'z': 0},
            'eye': {'x': -1.3, 'y': -1.6, 'z': 1.7}
        }
    )


    return figure


def visualize_mu_scan_E(df):

    figure = go.Figure()

    figure.update_layout(
        width=800, height=800,
        margin=dict(t=10, r=0, l=0, b=0)
    )


    # Add the (mu, R, S_z) scatter.
    R_values = df['R'].unique()
    for R in R_values:
        df_R = df[abs(df['R'] - R) < 1.0e-08]
        figure.add_trace(
            go.Scatter3d(
                x = df_R['mu'],
                y = df_R['R'],
                z = df_R['energy'],
                mode = 'markers+lines',
                marker = {
                    'size': 2,
                    'color': df_R['S_z'],
                    'colorscale': 'Sunsetdark'
                },
                line = {
                    'color': 'black',
                    'width': 0.2
                },
                showlegend = False
            )
        )


    # Add axis titles.
    figure.update_scenes(
        xaxis_title = "μ",
        yaxis_title = "Internuclear distance (a.u.)",
        zaxis_title = "E (a.u.)"
    )


    # Adjust the camera position.
    figure.update_layout(
        scene_camera = {
            'up': {'x': 0, 'y': 0, 'z': 1},
            'center': {'x': 0, 'y': 0, 'z': 0},
            'eye': {'x': -1.3, 'y': -1.6, 'z': 1.7}
        }
    )


    return figure
