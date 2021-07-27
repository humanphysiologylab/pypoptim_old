import pytest
import numpy as np
from ...algorythm.solution import Solution


@pytest.fixture()
def square_solution():

    class SquareSolution(Solution):
        def update(self):
            self._y = np.sum(self.x ** 2)

        def is_valid(self):
            return self.is_updated()

    return SquareSolution


class TestSolution:

    def test_init(self):

        invalid_xs = [42, 'foo', [], [[0, 1], [2, 3]]]
        for x in invalid_xs:
            with pytest.raises(ValueError):
                Solution(x)

        sol = Solution([42])
        assert not sol.is_updated()
        with pytest.raises(NotImplementedError):
            sol.is_valid()
        with pytest.raises(NotImplementedError):
            sol.update()

    def test_len(self):
        n = 42
        assert len(Solution(np.random.random(n))) == n

    def test_data(self):

        x = [0]

        sol = Solution(x)
        assert 'smth' not in sol
        with pytest.raises(KeyError):
            smth = sol['smth']

        kw = dict(a=1, b='foo')
        sol = Solution(x, **kw)
        for key, value in kw.items():
            assert sol[key] == value
            assert key in sol

        forty_two = 42
        sol['c'] = forty_two
        assert sol._data['c'] == forty_two

    def test_solution_update(self, square_solution):

        x = [1, 2]
        sol = square_solution(x)
        assert not sol.is_updated()
        assert not sol.is_valid()

        sol.update()
        assert sol.is_updated()
        assert sol.is_valid()
        assert sol.y == sum(xi**2 for xi in x)

        sol.x = x
        assert not sol.is_updated()
        assert not sol.is_valid()
