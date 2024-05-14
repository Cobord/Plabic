"""
test for double wiring diagrams
"""
#pylint:disable=invalid-name
from typing import List
from plabic import WiringDiagram

def test_dw1() -> None:
    """
    sbar_2 s_1 s2 sbar_1 sbar_2 s1
    """
    double_word = [-2,1,2,-1,-2,1]
    dw = WiringDiagram(3,double_word)
    assert dw.is_flip(True)
    assert dw.is_flip(False)
    expected_chamber_minors = [
        (set([1,2]),set([2,3])),
        (set([1]),set([3])),
        (set([1,2]),set([1,3])),
        (set([2]),set([3])),
        (set([2,3]),set([1,3])),
        (set([2]),set([1])),
        (set([3]),set([1])),
        (set([2,3]),set([1,2])),
        (set([1,2,3]),set([1,2,3]))
    ]
    for expected_chamber_minor,chamber_minor in zip(expected_chamber_minors,dw.chamber_minors()):
        assert chamber_minor == expected_chamber_minor
    p = dw.to_plabic()
    assert p.my_extra_props == set(["position"])
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert did_scale_up

def test_dw2() -> None:
    """
    empty word
    """
    DID_ERROR = False
    try:
        double_word : List[int] = []
        _dw = WiringDiagram(2,double_word)
    except ValueError:
        DID_ERROR = True
    assert DID_ERROR

def test_dw3() -> None:
    """
    s_1
    """
    double_word = [1]
    dw = WiringDiagram(2,double_word)
    assert dw.is_flip(True)
    assert not dw.is_flip(False)
    p = dw.to_plabic()
    assert p.is_bipartite()
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching == set([("wire1,0","wire2,0",1),("wire2,0","wire1,0",2)])
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert did_scale_up

def test_dw4() -> None:
    """
    s_1 sbar_1
    """
    double_word = [1,-1]
    dw = WiringDiagram(2,double_word)
    p = dw.to_plabic()
    assert dw.is_flip(True)
    assert dw.is_flip(False)
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching == set([("wire1,0","wire2,0",1),
                                         ("wire2,0","wire1,0",2),
                                         ("wire1,1","wire2,1",1),
                                         ("wire2,1","wire1,1",2)])
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*2,z[1]*2))
    assert did_scale_up
