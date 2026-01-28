# buildutils.py - utilities for build module


# _declass used to extract values from Property objects (legacy System class)
# With new AtomicScaleSample architecture, Property class is removed
# Keep _declass for now as a passthrough for backward compatibility
def _declass(item):
    """
    Pass through function for backward compatibility.
    Previously extracted value from Property objects.
    """
    return item
