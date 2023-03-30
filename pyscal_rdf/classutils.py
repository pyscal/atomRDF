class AttrSetter:
    def __init__(self):
        self._map_dict = {}
        self.head = None
    
    def __dir__(self):
        return list(self._map_dict.keys())

    def __repr__(self):
        return ", ".join(list(self._map_dict.keys()))    
    
    def __getattr__(self, key):
        """
        Calls attributes from the class when general lookup fails. There are four main types of
        calls that can occur.

        Normal attribute: this function does not get called
        key which is present in `_map_dict`
        """
        if key in self._map_dict.keys():
            if self.head is None:
                self.head = self
            return self._check_filter(key)
    
    def _check_filter(self, key):
        """
        Return the corresponding key
        """
        if isinstance(self._map_dict[key], str):
            return self.head[self._map_dict[key]]
        else:
            return self._map_dict[key]

    def _add_attribute(self, indict, head=None):
        """
        Add attributes the class in form of a dict
        """
        if head is None:
            head = self
        self.head = head

        for key, val in indict.items():
            if isinstance(val, dict):
                if key not in self._map_dict.keys():
                    self._map_dict[key] = AttrSetter()

                self._map_dict[key]._add_attribute(val, head=head)
            else:
                self._map_dict[key] = val