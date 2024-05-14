"""
test planar diagrams
"""
#pylint:disable=invalid-name,too-many-locals
from typing import cast,List,Tuple
from sympy import symbols, Symbol, Expr, Add, Mul, Pow, Integer, Rational, UnevaluatedExpr
from plabic import PlanarNetwork, determinant, ChipType

def test_planar() -> None:
    """
    example is taken from https://arxiv.org/pdf/math/9912128.pdf figure 1
    """

    def annotated_symbols(*args, **kwargs) -> Tuple[Symbol, ...]:
        """
        because symbols is typed Any
        despite the fact that for the arguments we are using it returns an tuple of Symbols
        """
        return symbols(*args, **kwargs)
    A, B, C, D, E, F, G, H, I = annotated_symbols(
        'a,b,c,d,e,f,g,h,i', commutative=True)
    ZERO, ONE = annotated_symbols('zero,one', commutative=True)
    p = PlanarNetwork(3, edge_list=[[(3, 2, A), (3, 3, ONE)], [(3, 2, C), (2, 1, B)], [
                      (2, 2, E), (1, 1, D)], [(3, 3, F), (2, 3, G), (1, 2, H)], [(2, 3, I)]],
                      multiplicative_identity=ONE,
                      additive_identity=ZERO)
    expected_weight_matrix: List[List[Expr]] = [
        [D, D*H, D*H*I], [B*D, B*D*H+E, B*D*H*I+E*(G+I)],\
            [A*B*D, A*B*D*H+(A+C)*E, A*B*D*H*I+(A+C)*E*(G+I)+F]]
    for i in range(1, 4):
        for j in range(1, 4):
            w_ij = cast(Expr, p.weight_matrix(i, j))
            w_ij = w_ij.subs({ZERO: 0.0, ONE: 1.0})
            assert w_ij.equals(expected_weight_matrix[i-1][j-1])
    expected_23_vertex_disjoint_weights = [B*C*D*E*G*H, E*F, B*D*F*H]
    for idx,weight in enumerate(p.vertex_disjoint_collection({2, 3}, {2, 3})):
        weight = cast(Expr, weight)
        weight = weight.subs({ZERO: 0.0, ONE: 1.0})
        assert weight.equals(expected_23_vertex_disjoint_weights[idx])
    delta_23_23 = cast(Expr, p.lindstrom_minor({2, 3}, {2, 3}))
    delta_23_23 = delta_23_23.subs({ZERO: 0.0, ONE: 1.0})
    expected = (B*C*D*E*G*H + B*D*F*H + E*F)*1.0
    assert delta_23_23 == expected
    slice_23 = (2, 3)
    my_minor = [[expected_weight_matrix[cur_i-1][cur_j-1]
                 for cur_i in slice_23] for cur_j in slice_23]
    expected_2 = determinant(my_minor, ONE,
                             lambda z: z.subs(
                                 {ZERO: 0.0, ONE: 1.0}).simplify(),
                             ZERO, -ONE)
    assert expected_2 == expected

    delta_123_123 = cast(Expr, p.lindstrom_minor({2, 3, 1}, {1, 2, 3}))
    delta_123_123 = delta_123_123.subs({ZERO: 0.0, ONE: 1.0})
    expected_determinant = determinant(expected_weight_matrix,
                                       ONE,
                                       lambda z: z.subs(
                                           {ZERO: 0.0, ONE: 1.0}).simplify(),
                                       ZERO,
                                       -ONE)
    assert delta_123_123 == expected_determinant
    def nn_with_my_symbols(edge_weight : Expr) -> bool:
        """
        a check that edge weight is nonnegative
        provided that the given symbols are nonnnegative
        checks with subtraction free
        """
        if isinstance(edge_weight,(Integer,Rational)):
            return edge_weight>0
        if isinstance(edge_weight,Symbol):
            return edge_weight in {ONE, A, B, C, D, E, F, G, H, I}
        op = edge_weight.func
        args = edge_weight.args
        if op == UnevaluatedExpr:
            return nn_with_my_symbols(cast(Expr,args[0]))
        if op in [Add,Mul]:
            return all(nn_with_my_symbols(cast(Expr,arg)) for arg in args)
        if op == Pow:
            return nn_with_my_symbols(cast(Expr,args[0])) and isinstance(args[1],(Integer,Rational))
        return False

    tnn = p.totally_nonnegative(nn_with_my_symbols)
    assert tnn
    assert p.chip_type == [(ChipType.DOWN,2),(ChipType.FLAT,3),\
                           (ChipType.DOWN,2),(ChipType.DOWN,1),\
                           (ChipType.FLAT,2),(ChipType.FLAT,1),\
                           (ChipType.FLAT,3),(ChipType.UP,3),\
                           (ChipType.UP,2),(ChipType.UP,3)]
