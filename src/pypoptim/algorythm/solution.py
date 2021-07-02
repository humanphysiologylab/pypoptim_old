class Solution:

    def __init__(self, x=None):
        self._x = x
        self._y = None
        self._data = dict()

    def __repr__(self):
        s =  'Solution = {\n'
        s += '    x    = {},\n'.format(self._x)
        s += '    y    = {},\n'.format(self._y)
        s += '    data = {}\n'.format(self._data)
        s += '}'
        return s

    def __str__(self):
        return self.__repr__()


    def __eq__(self, other):
        if not (self.is_valid() and other.is_valid()):
            raise ValueError("Can't compare invalid Solutions")
        return self.y == other.y

    def __gt__(self, other):
        if not (self.is_valid() and other.is_valid()):
            raise ValueError("Can't compare invalid Solutions")
        return self.y > other.y

    def __ne__(self, other):
        return not self.__eq__(other)

    def __le__(self, other):
        return not self.__gt__(other)

    def __ge__(self, other):
        return self.__gt__(other) or self.__eq__(other)

    def __lt__(self, other):
        return not self.__ge__(other)


    def __contains__(self, item):
        return item in self.data

    def __getitem__(self, item):
        if not self.__contains__(item):
            raise KeyError
        else:
            return self.data[item]

    def __setitem__(self, key, value):
        self._data[key] = value

    def __len__(self):
        return len(self._x)


    # x
    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self._x = x

    # y
    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        self._y = y

    # data
    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        self._data = data


    def update(self, *args, **kwargs) -> None:
        raise NotImplementedError("You must implement this method on your side!")

    def is_valid(self) -> bool:
        raise NotImplementedError("You must implement this method on your side!")
