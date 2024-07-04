"""
test for symbol alphabet
"""
from __future__ import annotations
from itertools import chain
from typing import Dict, List, Tuple

from plabic.cluster import Arithmetic
from plabic.ownership_workarounds import SameObjectError
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


def make_my_symbol_word(comma_separated: str) -> SymbolWord[TempExpression,TempExpression]:
    """
    make symbol word with tensorands being TempExpression
    """
    letters = comma_separated.split(",")
    return SymbolWord([TempExpression(letter) for letter in letters],
                      lambda z: z.vars(),
                      lambda z: z.dlog(),
                      TempExpression("a")**0)


def test_0log():
    """
    zero tensorands
    """
    errored = False
    try:
        _nonexistent = SymbolWord([],
                                  lambda z: z.vars(),
                                  lambda z: z.dlog(),
                                  TempExpression("a")**0)
    except ValueError:
        errored = True
    assert errored

def test_log():
    """
    symbol is just a single letter for a single log
    """
    loga = SymbolWord([TempExpression("a")],
                      lambda z: z.vars(),
                      lambda z: z.dlog(),
                      TempExpression("a")**0)
    assert loga.my_variables() == frozenset(["a"])
    da_over_a = loga.to_integrands()
    assert da_over_a == [{"a":TempExpression("a")**(-1)}]
    five_loga = Symbol([(5,loga),(2,loga)])
    assert len(five_loga) == 2

    logalpha = SymbolWord([TempExpression("alpha")],
                          lambda z: z.vars(),
                          lambda z: z.dlog(),
                          TempExpression("a")**0)
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
                           lambda z: z.dlog(),
                           TempExpression("a")**0)
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
                           lambda z: z.dlog(),
                           TempExpression("a")**0)
    assert logalphab1.am_i_zero()
    assert logalphab1.to_integrands() == [{"alpha":TempExpression("alpha")**(-1)},
                                          {"b":TempExpression("b")**(-1)},
                                          {"1":TempExpression("1")**(-1)}]
    zeros_out = Symbol([(4.3,logalphab1),(0,logalphab1),(2,logalphab1)])
    assert len(zeros_out) == 0

def test_disappeared_shuffle():
    """
    0*something
    and something*0
    """
    logalphab1 = SymbolWord([TempExpression("alpha"),TempExpression("b"),
                             TempExpression("dummy")**0],
                           lambda z: z.vars(),
                           lambda z: z.dlog(),
                           TempExpression("a")**0)
    assert logalphab1.am_i_zero()
    zeros_out = Symbol([(4.3,logalphab1),(0,logalphab1),(2,logalphab1)])
    assert len(zeros_out) == 0
    anything_else = SymbolWord([TempExpression("alpha"),TempExpression("b")],
                           lambda z: z.vars(),
                           lambda z: z.dlog(),
                           TempExpression("a")**0)
    zeros_out.ishuffle_product(Symbol([(1.0,anything_else)]))
    assert len(zeros_out) == 0

    logalphab1 = SymbolWord([TempExpression("alpha"),TempExpression("b"),
                             TempExpression("dummy")**0],
                           lambda z: z.vars(),
                           lambda z: z.dlog(),
                           TempExpression("a")**0)
    assert logalphab1.am_i_zero()
    zeros_out = Symbol([(4.3,logalphab1),(0,logalphab1),(2,logalphab1)])
    assert len(zeros_out) == 0
    anything_else = SymbolWord([TempExpression("alpha"),TempExpression("b")],
                           lambda z: z.vars(),
                           lambda z: z.dlog(),
                           TempExpression("a")**0)
    anything_else = Symbol([(1.0,anything_else)])
    assert len(anything_else) == 1
    anything_else.ishuffle_product(zeros_out)
    assert len(anything_else) == 0

def test_regular_shuffle():
    """
    factor_1 = 5 a otimes b + 3 b
    factor_2 = 10 x otimes y
    """
    #pylint:disable=unnecessary-lambda-assignment
    do_dlog = lambda z: z.dlog()
    do_vars = lambda z: z.vars()
    logab = SymbolWord([TempExpression("a"),TempExpression("b")],
                           do_vars,
                           do_dlog,
                           TempExpression("a")**0)
    logb = SymbolWord([TempExpression("b")],
                           do_vars,
                           do_dlog,
                           TempExpression("a")**0)
    logxy = SymbolWord([TempExpression("x"),TempExpression("y")],
                           do_vars,
                           do_dlog,
                           TempExpression("a")**0)
    factor_1 = Symbol([(5,logab),(3,logb),(0,logb)])
    assert len(factor_1) == 2
    factor_2 = Symbol([(10,logxy),(0,logab)])
    assert len(factor_2) == 1
    factor_1.ishuffle_product(factor_2)
    #pylint:disable=protected-access
    observed_answer = list(map(lambda summand: (summand[0],summand[1]._underlying),
                          list(factor_1._underlying)))
    expected_shuffles : List[Tuple[float,str]] = \
        list(map(lambda x: (50,x), ["a,b,x,y","a,x,b,y","x,a,b,y","a,x,y,b","x,a,y,b","x,y,a,b"]))
    expected_shuffles = \
        list(chain(expected_shuffles, map(lambda x: (30,x), ["x,y,b","x,b,y","b,x,y"])))
    for (expected_coeff, expected_summand) in expected_shuffles:
        expected_summand_word = make_my_symbol_word(expected_summand)
        assert (expected_coeff, expected_summand_word._underlying) in observed_answer
    assert len(factor_1) == len(expected_shuffles)

def test_clone():
    """
    repeat same symbol multiple times
    """
    logab = make_my_symbol_word("a,b")
    one_logab = Symbol([(1,logab),(0,logab)])
    assert len(one_logab) == 1
    repeated_list = [one_logab,one_logab,one_logab,one_logab]
    for i in range(4):
        for j in range(i):
            assert repeated_list[i] is repeated_list[j]
    #pylint:disable = unnecessary-lambda-assignment,protected-access
    temp_cloner = lambda z: TempExpression(str(z))
    nonrepeated_list = [z.clone(temp_cloner) for z in repeated_list]
    for i in range(4):
        for j in range(i):
            assert not nonrepeated_list[i] is nonrepeated_list[j]
            assert nonrepeated_list[i] != nonrepeated_list[j]
            are_same_i = list((both[0],both[1]._underlying)
                              for both in (nonrepeated_list[i])._underlying)
            are_same_j = list((both[0],both[1]._underlying)
                              for both in (nonrepeated_list[j])._underlying)
            assert are_same_i == are_same_j

def test_self_mul():
    """
    when trying to shuffle the same object with itself
    needs to clone one of them first
    if can't clone, then error
    if can then do so and give the expected answer
    """
    logab = make_my_symbol_word("a,b")
    one_logab = Symbol([(1,logab),(0,logab)])
    errored = False
    try:
        one_logab.ishuffle_product(one_logab)
    except SameObjectError as e:
        assert str(e) == "The objects must have distinct ids. We don't know how to clone for you."
        errored = True
    assert errored
    assert len(one_logab) == 1
    errored = False
    try:
        one_logab.ishuffle_many([one_logab])
    except SameObjectError as e:
        assert str(e) == "The objects must have distinct ids. We don't know how to clone for you."
        errored = True
    assert errored
    assert len(one_logab) == 1
    #pylint:disable=unnecessary-lambda-assignment
    temp_cloner = lambda z: TempExpression(str(z))

    # a a b b 4
    # a b a b 2
    expected_shuffles : List[Tuple[float,str]] = \
        [(2,"a,b,a,b"),(4,"a,a,b,b")]

    for_product = one_logab.clone(temp_cloner)
    for_product.ishuffle_product(for_product, how_to_clone_repeats=temp_cloner)
    for (expected_coeff, expected_summand) in expected_shuffles:
        expected_summand_word = make_my_symbol_word(expected_summand)
        #pylint:disable=protected-access
        coeff_this_word = sum(
            (coeff for (coeff, word) in for_product._underlying
             if word._underlying==expected_summand_word._underlying))
        assert coeff_this_word == expected_coeff
    assert len(for_product) == 6

    one_logab.ishuffle_many([one_logab], how_to_clone_repeats=temp_cloner)
    #pylint:disable=protected-access
    for (expected_coeff, expected_summand) in expected_shuffles:
        expected_summand_word = make_my_symbol_word(expected_summand)
        #pylint:disable=protected-access
        coeff_this_word = sum(
            (coeff for (coeff, word) in one_logab._underlying
             if word._underlying==expected_summand_word._underlying))
        assert coeff_this_word == expected_coeff
    assert len(for_product) == 6

def test_cell_var():
    """
    checking cell var from loop
    """
    #pylint:disable = import-outside-toplevel,cell-var-from-loop
    from functools import partial
    full_iterable = []
    use_partial = True
    for coeff in range(10):
        cur_coeff = coeff*4
        if use_partial:
            only_on_item = partial(
                lambda cur_coeff_2,summand: (cur_coeff_2,summand),
                cur_coeff)
            full_iterable = chain(full_iterable, map(only_on_item,range(4+coeff)))
        else:
            full_iterable = chain(full_iterable, map(lambda item: (cur_coeff,item),range(4+coeff)))
    assert list(full_iterable) == [(4*coeff,y) for coeff in range(10) for y in range(4+coeff)]
