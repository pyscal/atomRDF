import itertools
import numpy as np
from plotly import graph_objs as go

#color themes
from pyscal_rdf.gui.themes import themes

def create_box_plot(box, origin=[0,0,0]):
    """
    Create a plot of the simulation box
    """
    box = np.array(box)
    origin = np.array(origin)
    combos = list(itertools.combinations(range(3), 2))
    faces = []
    for combo in combos:
        f1 = [origin, box[combo[0]], box[combo[0]]+box[combo[1]], box[combo[1]], origin]
        s = combo[0] + combo[1]
        t = 3-s
        f2 = [origin + box[t], box[combo[0]]+ box[t],  box[combo[0]]+box[combo[1]]+ box[t], box[combo[1]]+ box[t], origin + box[t]]
        faces.append(np.array(f1))
        faces.append(np.array(f2))
    traces = []
    for face in faces:
        trace = go.Scatter3d(
            x=face[:,0],
            y=face[:,1],
            z=face[:,2],
            mode='lines',
            name='lines',
            line=dict(width=2.0, color='#263238'),
            showlegend=False
        )
        traces.append(trace)
    return traces

def plot_3d(sys, theme="teal"):
    """
    Plot 3D of simulation box
    """
    pos = np.array(sys.atoms.positions)
    traces = create_box_plot(sys.box)
    radius = 10
    opacity = 1.0
    color = themes[theme]["atom"]
    data=go.Scatter3d(
        x=pos[:,0],
        y=pos[:,1],
        z=pos[:,2],
        mode='markers',
        opacity=1.0,
        marker=dict(
            sizemode='diameter',
            sizeref=750,
            size=radius,
            color = color,
            opacity = opacity,
            line=dict(width=0.5, color='#455A64')
        ),
    )
    traces.append(data)
    fig = go.Figure(data=traces)
    fig.update_layout(scene = dict(
                        xaxis_title="",
                        yaxis_title="",
                        zaxis_title="",
                        xaxis = dict(
                             showticklabels=False,
                             showbackground=False,
                             zerolinecolor="#455A64",),
                        yaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64"),
                        zaxis = dict(
                            showticklabels=False,
                            showbackground=False,
                            zerolinecolor="#455A64",),),
                        width=300,
                        height=300,
                        margin=dict(
                        r=10, l=10,
                        b=10, t=10)
                      )
    fig.update_layout(showlegend=False)
    fig.show()


def system(sys, output, theme="teal"):
    output.clear_output()
    with output:
        plot_3d(sys, theme=theme)