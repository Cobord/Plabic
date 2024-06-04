"""
test for symbol alphabet
"""
from __future__ import annotations
from typing import Dict, List

from plabic.cluster import Arithmetic
from plabic.symbol_alphabet import Symbol,SymbolWord

class TempExpression(Arithmetic):
    """
    a string for a rational expression
    need a wrapper class so that it implements the
    Arithmetic protocol
    just enough of the implementation for the below tests
    would not do to_integrands properly on others
    """
    def __init__(self, me: str):
        self._unparsed = me

    def __add__(self,other : TempExpression) -> TempExpression:
        return TempExpression(f"({self}) + ({other})")

    def __mul__(self,other : TempExpression) -> TempExpression:
        return TempExpression(f"({self}) * ({other})")

    def __pow__(self,power : int) -> TempExpression:
        if power == 0:
            return TempExpression("1")
        return TempExpression(f"({self}) ^ ({power})")

    def __truediv__(self,other : TempExpression) -> TempExpression:
        return TempExpression(f"({self}) / ({other})")

    def __str__(self) -> str:
        return self._unparsed

    def __eq__(self, other : TempExpression) -> bool:
        return self._unparsed == other._unparsed

    def vars(self) -> List[str]:
        """
        list the variables contained in the expression
        """
        return [self._unparsed]

    def dlog(self) -> Dict[str,TempExpression]:
        """
        assumes that self is a single variable
        as in the below tests
        """
        return {str(self) : self**(-1)}



def test_0log():
    """
    zero tensorands
    """
    errored = False
    try:
        _nonexistent = SymbolWord([],
                                  lambda z: z.vars(),
                                  lambda z: z.dlog())
    except ValueError:
        errored = True
    assert errored

def test_log():
    """
    symbol is just a single letter for a single log
    """
    loga = SymbolWord([TempExpression("a")],
                      lambda z: z.vars(),
                      lambda z: z.dlog())
    assert loga.my_variables() == frozenset(["a"])
    da_over_a = loga.to_integrands()
    assert da_over_a == [{"a":TempExpression("a")**(-1)}]
    five_loga = Symbol([(5,loga),(2,loga)])
    assert len(five_loga) == 2

    logalpha = SymbolWord([TempExpression("alpha")],
                          lambda z: z.vars(),
                          lambda z: z.dlog())
    assert logalpha.my_variables() == frozenset(["alpha"])
    dalpha_over_alpha = logalpha.to_integrands()
    assert dalpha_over_alpha == [{"alpha":TempExpression("alpha")**(-1)}]
    eight_logalpha = Symbol([(8,logalpha),(0,logalpha),(0,logalpha),(0,logalpha)])
    assert len(eight_logalpha) == 1


def test_2log():
    """
    two tensorands
    int dalpha/alpha int db/b
    log alpha log b + C log alpha + D log b + CD
    """
    logalphab = SymbolWord([TempExpression("alpha"),TempExpression("b")],
                           lambda z: z.vars(),
                           lambda z: z.dlog())
    assert logalphab.my_variables() == frozenset(["alpha","b"])
    diff_forms = logalphab.to_integrands()
    assert diff_forms == [{"alpha":TempExpression("alpha")**(-1)},{"b":TempExpression("b")**(-1)}]

def test_1_disappear():
    """
    alpha otimes b otimes 1
    disappears
    """
    logalphab1 = SymbolWord([TempExpression("alpha"),TempExpression("b"),
                             TempExpression("dummy")**0],
                           lambda z: z.vars(),
                           lambda z: z.dlog())
    assert logalphab1.am_i_zero()
    zeros_out = Symbol([(4.3,logalphab1),(0,logalphab1),(2,logalphab1)])
    assert len(zeros_out) == 0

def test_cell_var():
    """
    checking cell var from loop
    """
    #pylint:disable = import-outside-toplevel
    import functools
    from itertools import chain
    full_iterable = []
    use_partial = True
    for coeff in range(10):
        cur_coeff = coeff*4
        if use_partial:
            only_on_item = functools.partial(
                lambda cur_coeff_2,summand: (cur_coeff_2,summand),
                cur_coeff)
            full_iterable = chain(full_iterable, map(only_on_item,range(4+coeff)))
        else:
            #pylint:disable=cell-var-from-loop
            full_iterable = chain(full_iterable, map(lambda item: (cur_coeff,item),range(4)))
    assert list(full_iterable) == [(4*coeff,y) for coeff in range(10) for y in range(4+coeff)]
