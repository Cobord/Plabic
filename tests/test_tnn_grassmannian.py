"""
test for tnn grassmannian
"""
from sympy import Integer,symbols
from plabic.tnn_grassmannian import TNNGrassChart,MinorSignConstraints

def test_gr24():
    """
    R^3_gt 0 subsets of Gr_2,4
    """
    #pylint:disable = invalid-name
    A,B,C = symbols("a,b,c")
    ZERO = Integer(0)
    ONE = Integer(1)
    positivity_constraints = MinorSignConstraints(
        [set([0,1]),set([0,2]),set([1,2]),set([1,3]),set([2,3])],
        [],[set([0,3])],
        [],[],[])
    gr24 = TNNGrassChart(my_matrix=[[ONE,ZERO,-A,-B],[ZERO,ONE,C,ZERO]],
                         variables=[A,B,C],
                         sign_constraints=positivity_constraints,
                         check_tnn=True)
    assert gr24.sign_constraints.is_positroidy(4)
    gr24.validate_pluckers()
    assert gr24.includes_in_which_schubert() == set([0,1])
    positivity_constraints = MinorSignConstraints(
        [set([0,2]),set([1,2]),set([1,3]),set([2,3]),set([0,3])],
        [],[set([0,1])],
        [],[],[])
    gr24 = TNNGrassChart(my_matrix=[[ONE,A,ZERO,-B],[ZERO,ZERO,ONE,C]],
                         variables=[A,B,C],
                         sign_constraints=positivity_constraints,
                         check_tnn=True,gives_positroid_cell=True)
    assert gr24.includes_in_which_schubert() == set([0,2])
    assert gr24.sign_constraints.is_positroidy(4)
    gr24.cyclic_shift(800)
    assert gr24.includes_in_which_schubert() == set([0,2])
    gr24.cyclic_shift(2001)
    assert gr24.includes_in_which_schubert() == set([0,1])


def test_gr25():
    """
    R^6_gt 0 subset of Gr_2,5
    """
    #pylint:disable = invalid-name
    A1,A2,A3,A4,A5,A6 = symbols("a1,a2,a3,a4,a5,a6")
    ZERO = Integer(0)
    ONE = Integer(1)
    positivity_constraints = MinorSignConstraints.all_positive(5,2)
    gr25 = TNNGrassChart(my_matrix=[[ONE,A1+A2*A6,A6,A3*A6,ZERO],
                                    [ZERO,A2*A5*A6,A5*A6,A4+A3*A5*A6,ONE]],
                         variables=[A1,A2,A3,A4,A5,A6],
                         sign_constraints=positivity_constraints,
                         check_tnn=True,gives_positroid_cell=True)
    assert gr25.sign_constraints.is_positroidy(5)
    assert gr25.includes_in_which_schubert() == set([0,1])
