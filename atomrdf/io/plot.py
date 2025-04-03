import pyscal3.operations.visualize as visualize
from atomrdf.io.write import write

def plot(system, plot_style='all', colorby=None, 
    cmap = 'viridis', 
    radius=10, 
    opacity=1.0,
    hide_zero=False,
    color = '#ff7f00'):
    
    if plot_style == 'all':
        visualize.plot_simple(system)
    
    elif plot_style == 'selection':
        visualize.plot_by_selection(system, radius=radius, opacity=opacity)
    
    elif plot_style == 'continuous_property':
        visualize.plot_by_property(system, colorby, 
            cmap = cmap, 
            radius=radius, 
            opacity=opacity,
            hide_zero=hide_zero)
        
    elif plot_style == 'boolean_property':
        visualize.plot_by_boolean(system, colorby, 
            color = color, 
            radius=radius, 
            opacity=opacity,
            hide_zero=hide_zero)

def plot3d(system, *args, **kwargs):
    try:
        from pyiron_atomistics.atomistics.structure.atoms import (
            ase_to_pyiron,
            pyiron_to_ase,
        )
    except ImportError:
        raise ImportError("Please install pyiron_atomistics")
    ase_structure = write(system, format='ase')
    pyiron_structure = ase_to_pyiron(ase_structure)
    return pyiron_structure.plot3d(*args, **kwargs)