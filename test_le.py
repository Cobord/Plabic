"""
test Le diagrams
"""

#pylint:disable=invalid-name,import-error
from typing import Set
from le_diagram import LeDiagram

def test_empty() -> None:
    """
    partition should be of positive number
    """
    try:
        _my_Le = LeDiagram([])
    except ValueError as e:
        assert str(e)=="Expect the shape to be a partition of a positive number"
    try:
        _my_Le = LeDiagram([[],[],[]])
    except ValueError as e:
        assert str(e)=="Expect the shape to be a partition of a positive number"

def test_Le1() -> None:
    """
    Example from
    https://doi.org/10.1016/j.jcta.2018.06.001
    """
    my_diagram = [[0,1,0,1,0],[1,1,0,1],[0,0],[0,1]]
    my_Le = LeDiagram(my_diagram)
    assert my_Le.filling == my_diagram
    print(my_Le.column_height(0) == 4)
    print(my_Le.column_height(1) == 4)
    print(my_Le.column_height(2) == 2)
    print(my_Le.column_height(3) == 2)
    print(my_Le.column_height(4) == 1)
    print(my_Le.column_height(5) == 0)
    p = my_Le.to_plabic()
    assert p.my_extra_props == set(["position"])
    assert p.my_perfect_matching is None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert did_scale_up

def test_rectangle() -> None:
    """
    case of n by k rectangle with all 1's
    have particularly nice face quivers
    """
    k,n = 5,4
    my_diagram = [[1]*k]*n
    my_Le = LeDiagram(my_diagram)
    p = my_Le.to_plabic()
    assert p.my_extra_props == set(["position"])
    assert p.my_perfect_matching is None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert did_scale_up

def from_multidigit_int(multi_digit : int) -> Set[int]:
    """
    when n is <=9
    can give the subset of 1..n
    as just concatenation
    """
    return {int(i) for i in str(multi_digit)}

def test_to_Grassmann() -> None:
    """
    example to Grassmann necklace
    from https://arxiv.org/pdf/1803.01726.pdf
    """
    my_diagram = [[1,1,0,0,1],[1,1,0,1],[0,1]]
    my_Le = LeDiagram(my_diagram)
    obs_necklace = my_Le.to_grassmann_necklace(3,8)
    exp_necklace_nums = [136,236,367,467,678,678,178,168]
    exp_necklace = [from_multidigit_int(z) for z in exp_necklace_nums]
    exp_necklace_2 = [set([1,3,6]),set([2,3,6]),set([3,6,7]),
                    set([4,6,7]),set([6,7,8]),set([6,7,8]),
                    set([1,7,8]),set([1,6,8])]
    assert exp_necklace == exp_necklace_2
    assert obs_necklace == exp_necklace
    p = my_Le.to_plabic()
    assert p.my_extra_props == set(["position"])
    assert p.my_perfect_matching is None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert did_scale_up

def test_from_Grassmann() -> None:
    """
    example from Grassmann necklace
    from https://arxiv.org/pdf/1803.01726.pdf
    """
    necklace_nums = [1247, 2347, 3478, 4678, 5678, 4678, 1478, 1478]
    necklace = [from_multidigit_int(z) for z in necklace_nums]
    print(necklace)
    my_Le = LeDiagram.from_grassmann_necklace(necklace,4,8)
    print(my_Le)
    exp_Le = LeDiagram([[1,0,0,1],[1,1,0,1],[0,0,1],[0]])
    assert my_Le.filling == exp_Le.filling
    p = my_Le.to_plabic()
    assert p.my_extra_props == set(["position"])
    assert p.my_perfect_matching is None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert did_scale_up
