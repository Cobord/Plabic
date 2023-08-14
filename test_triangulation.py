"""
test triangulations of convex m-gons
"""
#pylint:disable=import-error,invalid-name,R0801
from math import pi as PI,sin,cos
from triangulation import Triangulation

def test_octagon() -> None:
    """
    a triangulation of an octagon
    """
    N = 8
    my_diagonals = [(4,7),(2,8),(2,4),(2,7),(4,6)]
    RADIUS = 2
    t = Triangulation([ (RADIUS*cos(2*PI*which/N),RADIUS*sin(2*PI*which/N))
                       for which in range(N)], [(x-1,y-1) for x,y in my_diagonals])
    diags_before = t.my_diagonal_list.copy()
    for diag,quad in zip(t.my_diagonal_list,t.my_quadrilaterals):
        assert diag[0]==quad[0]
        assert diag[1]==quad[2]
    p = t.to_plabic()
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.10,z[1]*1.1))
    assert did_scale_up
    change, new_diag = t.quad_flip((1,3))
    assert change
    assert new_diag == (2,6)
    for which_diag_now,diag in enumerate(t.my_diagonal_list):
        if diag != (2,6):
            assert diag == diags_before[which_diag_now]
        else:
            assert diags_before[which_diag_now] == (1,3)
    for diag,quad in zip(t.my_diagonal_list,t.my_quadrilaterals):
        assert diag[0]==quad[0]
        assert diag[1]==quad[2]
    p = t.to_plabic()
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.3,z[1]*2.4))
    assert did_scale_up

def test_pentagon() -> None:
    """
    a triangulation of an pentagon
    """
    N = 5
    my_diagonals = [(1,3),(1,4)]
    RADIUS = 2
    t = Triangulation([ (RADIUS*cos(2*PI*which/N),RADIUS*sin(2*PI*which/N))
                       for which in range(N)], [(x-1,y-1) for x,y in my_diagonals])
    diags_before = t.my_diagonal_list.copy()
    for diag,quad in zip(t.my_diagonal_list,t.my_quadrilaterals):
        assert diag[0]==quad[0]
        assert diag[1]==quad[2]
    p = t.to_plabic()
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.50,z[1]*1.1))
    assert did_scale_up
    change, new_diag = t.quad_flip((0,2))
    assert change
    assert new_diag == (1,3)
    for which_diag_now,diag in enumerate(t.my_diagonal_list):
        if diag != (1,3):
            assert diag == diags_before[which_diag_now]
        else:
            assert diags_before[which_diag_now] == (0,2)
    for diag,quad in zip(t.my_diagonal_list,t.my_quadrilaterals):
        assert diag[0]==quad[0]
        assert diag[1]==quad[2]
    p = t.to_plabic()
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.3,z[1]*2.4))
    assert did_scale_up
