"""
Wrapper around Materials Project to query structures and get it as a KG
"""

from mp_api.client import MPRester
import numpy as np

def query_mp(api_key, chemical_system=None, material_ids=None, is_stable=True):
    rest = {
            "use_document_model": False,
            "include_user_agent": True,
            "api_key": api_key, 
        }
    if (chemical_system is None) and (material_ids is None):
        raise ValueError("Please provide either a chemical system or a list of material ids")
    with MPRester(**rest) as mpr:
        if chemical_system is not None:
            docs = mpr.materials.summary.search(chemsys=chemical_system, is_stable=is_stable)
        else:
            docs = mpr.materials.summary.search(material_ids=material_ids)
    return docs