class SortableObject(object):
    """Base class for GillesPy2 objects that are sortable."""

    def __eq__(self, other):
        return str(self)==str(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def __lt__(self, other):
        if hasattr(self, 'id') and hasattr(other, 'id'):
            return self.id.lower() < other.id.lower()
        elif hasattr(self, 'name') and hasattr(other, 'name'):
            return self.name.lower() < other.name.lower()
        else:
            return repr(self) < repr(other)

    def __le__(self, other):
        if hasattr(self, 'id') and hasattr(other, 'id'):
            return self.id.lower() <= other.id.lower()
        elif hasattr(self, 'name') and hasattr(other, 'name'):
            return self.name.lower() <= other.name.lower()
        else:
            return repr(self) <= repr(other)

    def __cmp__(self, other):
        if hasattr(self, 'id') and hasattr(other, 'id'):
            return cmp(self.id.lower(), other.id.lower())
        elif hasattr(self, 'name') and hasattr(other, 'name'):
            return cmp(self.name.lower(), other.name.lower())
        else:
            return cmp(repr(self), repr(other))

    def __hash__(self):
        if hasattr(self, '_hash'):
            return self._hash
        if hasattr(self, 'id'):
            self._hash = hash(self.id)
        elif hasattr(self, 'name'):
            self._hash = hash(self.name)
        else:
            self._hash = hash(self)
        return self._hash