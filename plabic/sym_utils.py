"""
utilities with symbolic expressions
"""

from typing import cast,Set,List,Optional
from sympy import Symbol, Expr, Add, Mul, Pow, Integer, Rational, UnevaluatedExpr

#pylint:disable=too-many-return-statements
def nn_with_my_symbols(edge_weight : Expr,
                       must_be_positive : bool,
                       positive_symbols : Set[Symbol]) -> bool:
    """
    a check that edge weight is nonnegative
    provided that the given symbols are nonnnegative
    checks with subtraction free
    """
    if isinstance(edge_weight,(Integer,Rational)):
        if must_be_positive:
            return edge_weight>0
        return edge_weight>=0
    if isinstance(edge_weight,Symbol):
        return edge_weight in positive_symbols
    sym_op = edge_weight.func
    args = edge_weight.args
    if sym_op == UnevaluatedExpr:
        return nn_with_my_symbols(cast(Expr,args[0]),must_be_positive,positive_symbols)
    if sym_op == Add:
        all_nn = all(nn_with_my_symbols(cast(Expr,arg),False,positive_symbols) for arg in args)
        if must_be_positive:
            extra_condition = \
                any(nn_with_my_symbols(cast(Expr,arg),True,positive_symbols) for arg in args)
        else:
            extra_condition = True
        return all_nn and extra_condition
    if sym_op == Mul:
        all_nn = all(nn_with_my_symbols(cast(Expr,arg),False,positive_symbols) for arg in args)
        if must_be_positive:
            extra_condition = not any(arg==0 for arg in args)
        else:
            extra_condition = True
        return all_nn and extra_condition
    if sym_op == Pow:
        if args[1]==0:
            return True
        base_nonnegativity = nn_with_my_symbols(cast(Expr,args[0]),
                                                must_be_positive,positive_symbols)
        return base_nonnegativity and \
            isinstance(args[1],(Integer,Rational))
    return False

def determinant(matrix : List[List[Expr]], mul : Expr = Integer(1),
                simplify_freq : int=3,
                original_width : Optional[int] =None) -> Expr:
    """
    determinant using only + and *
    so that it can work on complex object types
    that override __add__ and __mul__
    imported det functions assume primitive numeric type like float
    """
    width = len(matrix)
    if original_width is None:
        original_width = width
    if width == 1:
        to_return = cast(Expr,mul * matrix[0][0])
        if original_width == width:
            to_return = cast(Expr,to_return.simplify())
        return to_return
    sign : Integer = Integer(-1)
    answer = cast(Expr,Integer(0))
    for skip_col in range(width):
        temp_matrix : List[List[Expr]] = []
        for jdx in range(1, width):
            buff : List[Expr] = []
            for kdx in range(width):
                if kdx != skip_col:
                    buff.append(matrix[jdx][kdx])
            temp_matrix.append(buff)
        sign *= Integer(-1)
        cur_det = determinant(temp_matrix, sign * matrix[0][skip_col],
                              simplify_freq, original_width)
        answer = cast(Expr,answer + mul * cur_det)
    if width % simplify_freq == original_width % simplify_freq:
        answer = cast(Expr,answer.simplify())
    return answer
