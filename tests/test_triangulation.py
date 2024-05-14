"""
test triangulations of convex m-gons
"""
#pylint:disable=invalid-name,R0801
from math import pi as PI,sin,cos
from plabic import Triangulation

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
    scale_factor = 2
    assert p.circles_config is not None
    assert len(p.circles_config.internal_circles)==0
    center, radius, _ = p.circles_config.outer_circle
    assert center[0]**2+center[1]**2 <= 1e-6
    assert radius == RADIUS*scale_factor
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.10,z[1]*1.1))
    assert did_scale_up
    assert p.circles_config is None
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
    scale_factor = 2
    assert p.circles_config is not None
    assert len(p.circles_config.internal_circles)==0
    center, radius, _ = p.circles_config.outer_circle
    assert center[0]**2+center[1]**2 <= 1e-6
    assert radius == RADIUS*scale_factor
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.3,z[1]*2.4))
    assert did_scale_up
    scale_factor = 2
    assert p.circles_config is None

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
    scale_factor = 2
    assert p.circles_config is not None
    assert len(p.circles_config.internal_circles)==0
    center, radius, _ = p.circles_config.outer_circle
    assert center[0]**2+center[1]**2 <= 1e-6
    assert radius == RADIUS*scale_factor
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.50,z[1]*1.1))
    assert did_scale_up
    assert p.circles_config is None
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
    scale_factor = 2
    assert p.circles_config is not None
    assert len(p.circles_config.internal_circles)==0
    center, radius, _ = p.circles_config.outer_circle
    assert center[0]**2+center[1]**2 <= 1e-6
    assert radius == RADIUS*scale_factor
    assert p.my_extra_props == set(["position","my_perfect_edge"])
    assert p.my_perfect_matching is not None
    did_scale_up = p.coordinate_transform(lambda z: (z[0]*1.3,z[1]*2.4))
    assert did_scale_up
    assert p.circles_config is None
