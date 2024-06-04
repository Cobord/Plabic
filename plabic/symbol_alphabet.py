"""
A symbol alphabet is TODO
"""
from __future__ import annotations
import functools
from itertools import chain
from math import comb
from multiprocessing import Pool
from typing import Callable, Dict, FrozenSet, Generic, Iterable, List, Set, Tuple, TypeVar

from .cluster import Arithmetic
from .shuffle import shuffle_product

Tensorand = TypeVar("Tensorand", bound=Arithmetic)
T = TypeVar("T")

#pylint:disable=too-few-public-methods
class SymbolWord(Generic[Tensorand,T]):
    """
    a_1 otimes ... a_m
    where a_i are all of type Tensorand
    there is a provided function which gives the variables in each of
        these a_i
    dlog tells how to turn each tensorand into da_i/a_i
        given as a dictionary from the variable names to the
        coefficient of d{that var}
        the type of such expressions is T which
        will often be the same as type as Tensorand
        both being Expression Trees for rational functions in those variables
    """
    def __init__(self,
                 letters: List[Tensorand],
                 letters_variables: Callable[[Tensorand],Iterable[str]],
                 dlog : Callable[[Tensorand],Dict[str,T]]
                 ):
        if len(letters) == 0:
            raise ValueError("There should be at least one tensorand")
        self._underlying = letters
        self.__used_letters : Set[str] = set()
        self.__no_variable_idces = []
        for (cur_idx,letter) in enumerate(letters):
            vars_in_this_letter = 0
            for var in letters_variables(letter):
                vars_in_this_letter = 0
                self.__used_letters.add(var)
            if vars_in_this_letter == 0:
                self.__no_variable_idces.append(cur_idx)
        self._dlog = dlog
        self._letters_variables = letters_variables

    def to_integrands(self) -> List[Dict[str,T]]:
        """
        da/a, db/b, ... each one as an element of type T
        """
        return [self._dlog(z) for z in self._underlying]

    def my_variables(self) -> FrozenSet[str]:
        """
        what variables are used
        """
        return frozenset(self.__used_letters)

    def am_i_zero(self) -> bool:
        """
        one of the tensorands is the multiplicative identity
        """
        one_tensorand = self._underlying[0]**0
        return one_tensorand in self._underlying

    def word_shuffle_product(self, other: SymbolWord):
        """
        summands in shuffle product
        """
        #pylint:disable=protected-access
        if self._dlog != other._dlog:
            raise ValueError("They should have the same way of translating tensorands f to dlog(f)")
        if self._letters_variables != other._letters_variables:
            raise ValueError("They should have the same way of translating tensorands f to dlog(f)")
        for z in shuffle_product(self._underlying, other._underlying):
            yield SymbolWord(z, self._letters_variables, self._dlog)

    def count_word_shuffle_product(self, other: SymbolWord) -> int:
        """
        summands in shuffle product
        """
        n1 = len(self._underlying)
        #pylint:disable=protected-access
        n2 = len(other._underlying)
        return comb(n1+n2,n1)

class Symbol:
    """
    a linear combination of SymbolWord
    the underlying summands only has to be an iterable
    with all the summands not a list
    """
    def __init__(self, summands: List[Tuple[float,SymbolWord]]):
        self._underlying : Iterable[Tuple[float,SymbolWord]] = [(coeff,summand) for (coeff,summand)
                             in summands if (coeff != 0 and not summand.am_i_zero())]
        self._num_summands = len(self._underlying)

    def __len__(self) -> int:
        return self._num_summands

    def ishuffle_many(self, rest: List[Symbol],*,
                      materialize_length = 200,
                      do_parallel=True,
                      how_many_more_times_do_parallel = 2):
        """
        replace self with shuffle product of self and all the rest
        using the associativity in order to split the computation
        if the number of factors is too large
        all of the original Symbol's are consumed in the process
        TODO test
        """
        num_factors = len(rest)+1
        if num_factors<5 or not do_parallel:
            for other in rest:
                self.ishuffle_product(other, materialize_length=materialize_length)
        else:
            split_half = num_factors//2
            with_self = rest[0:split_half]
            without_self_head = rest[split_half]
            without_self_tail = rest[(split_half+1):]
            if how_many_more_times_do_parallel > 0:
                do_parallel_next = True
                how_many_more_times_do_parallel_next = how_many_more_times_do_parallel-1
            else:
                do_parallel_next = False
                how_many_more_times_do_parallel_next = 0
            with Pool(2) as pool:
                both_results = [pool.apply_async(lambda main,remaining:
                                 main.ishuffle_many(remaining,
                                                    materialize_length=materialize_length,
                                                    do_parallel = do_parallel_next,
                                                    how_many_more_times_do_parallel=\
                                                        how_many_more_times_do_parallel_next
                                                    ),
                                 cur_arg) for cur_arg in
                                 [(self,with_self),
                                  (without_self_head,without_self_tail)]]
                _ = [result.get() for result in both_results]
            self.ishuffle_product(without_self_head,materialize_length=materialize_length)

    def ishuffle_product(self, other: Symbol, *, materialize_length: int = 200):
        """
        replaces self with shuffle product of self and other
        only materializes into a list instead of lazy generator chaining
        if the length is small enough
        other is consumed in the process
        TODO test
        """
        #pylint:disable=protected-access
        self_summands = []
        self_summands, self._underlying = self._underlying, self_summands
        self._num_summands = 0
        for (coeff_self,word_self) in self_summands:
            for (coeff_other,word_other) in other._underlying:
                cur_coeff = coeff_self*coeff_other
                if cur_coeff == 0:
                    # know from construction that coeff_self and coeff_other are nonzero
                    # so this shouldn't matter
                    continue
                num_shuffles = word_self.count_word_shuffle_product(word_other)
                self._underlying = chain(self._underlying,
                                         map(
                                             functools.partial(
                                                 lambda cur_coeff,summand: (cur_coeff,summand),
                                                 cur_coeff),
                                             word_self.word_shuffle_product(word_other)))
                self._num_summands += num_shuffles
                # know from construction that word_self and word_other am_i_zero()
                # are both giving false
                # so the summands should also have am_i_zero give false
                # could do the check with a continue like with cur_coeff, but didn't here
        if self._num_summands<materialize_length:
            self._underlying = list(self._underlying)
