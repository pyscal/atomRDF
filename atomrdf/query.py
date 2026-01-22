import numpy as np
import warnings
from atomrdf.build.bulk import _generate_atomic_sample_data
from atomrdf.datamodels.structure import AtomicScaleSample


def materials_project(
    api_key,
    chemical_system=None,
    material_ids=None,
    is_stable=True,
    conventional=True,
    graph=None,
):
    """
    Fetch structure(s) from Materials Project database.

    .. deprecated::
        Use :meth:`AtomicScaleSample.from_repository` instead.
        This function will be removed in a future version.

    Parameters
    ----------
    api_key : str
        Materials Project API key
    chemical_system : str, optional
        Chemical system string (e.g., 'Fe-C')
    material_ids : list of str, optional
        List of mp-ids (e.g., ['mp-149'])
    is_stable : bool, optional
        Only fetch stable materials (default: True)
    conventional : bool, optional
        Use conventional cell (default: True)
    graph : KnowledgeGraph, optional
        Knowledge graph to add structures to

    Returns
    -------
    AtomicScaleSample or list
        Single sample or list of samples

    See Also
    --------
    AtomicScaleSample.from_repository : Recommended replacement method
    """
    warnings.warn(
        "materials_project() is deprecated. Use AtomicScaleSample.from_repository() instead.",
        DeprecationWarning,
        stacklevel=2,
    )

    return AtomicScaleSample.from_repository(
        repository="materials_project",
        api_key=api_key,
        material_ids=material_ids,
        chemical_system=chemical_system,
        is_stable=is_stable,
        conventional=conventional,
        graph=graph,
    )
