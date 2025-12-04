from atomrdf.sample import Property

#declassing special variables
def _declass(item):
    if isinstance(item, Property):
        return item.value
    else:
        return item
