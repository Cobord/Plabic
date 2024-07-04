"""
symbol alphabet
as in the context of amplitudes
"""
from __future__ import annotations
import functools
from itertools import chain
from math import comb
from multiprocessing import Pool
from typing import Callable, Dict, FrozenSet, Generic, Iterable, List, Optional, Set, Tuple, TypeVar

from .ownership_workarounds import ConsumedObjectError, SameObjectError,\
    avoid_common_pointers, avoid_this_pointer

from .shuffle import shuffle_product

Tensorand = TypeVar("Tensorand")
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
    T will often be the same as type as Tensorand
        both being Expression Trees for rational functions in those variables
    """
    def __init__(self,
                 letters: List[Tensorand],
                 letters_variables: Callable[[Tensorand],Iterable[str]],
                 dlog : Callable[[Tensorand],Dict[str,T]],
                 makes_word_zero: Optional[Tensorand] = None
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
        self.__makes_word_zero = makes_word_zero

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
        if self.__makes_word_zero is not None:
            return self.__makes_word_zero in self._underlying
        return False

    def word_shuffle_product(self, other: SymbolWord):
        """
        summands in shuffle product
        """
        #pylint:disable=protected-access,pointless-string-statement
        if self._dlog != other._dlog:
            """
            should we assume that this is just two functions that
            do the same thing, but not literally the same?
            """
            raise ValueError("They should have the same way of translating tensorands f to dlog(f)")
        if self._letters_variables != other._letters_variables:
            """
            should we assume that this is just two functions that
            do the same thing, but not literally the same?
            """
            raise ValueError(
                "They should have the same way of saying which variables tensorands contain")
        for z in shuffle_product(self._underlying, other._underlying):
            yield SymbolWord(z, self._letters_variables, self._dlog, self.__makes_word_zero)

    def count_word_shuffle_product(self, other: SymbolWord) -> int:
        """
        summands in shuffle product
        """
        n1 = len(self._underlying)
        #pylint:disable=protected-access
        n2 = len(other._underlying)
        return comb(n1+n2,n1)

    def _def_nonzero(self) -> SymbolWord:
        """
        when already checked the summand is nonzero
        this extra property doesn't make a difference anymore
        """
        self.__makes_word_zero = None
        return self

    def clone(self, letter_copy: Callable[[Tensorand],Tensorand]) -> SymbolWord:
        """
        clone
        """
        letters_copy = [letter_copy(letter) for letter in self._underlying]
        return SymbolWord(letters_copy,
                               self._letters_variables,
                               self._dlog,
                               self.__makes_word_zero)

class Symbol(Generic[Tensorand,T]):
    """
    a linear combination of SymbolWord
    the underlying summands only has to be an iterable
    with all the summands not a list
    """
    def __init__(self, summands: List[Tuple[float,SymbolWord[Tensorand,T]]]):
        self._underlying : Iterable[Tuple[float,SymbolWord[Tensorand,T]]] = \
            [(coeff,summand._def_nonzero())
             for (coeff,summand) in summands
             if (coeff != 0 and not summand.am_i_zero())]
        self._num_summands = len(self._underlying)
        self._invalidated = False

    def __len__(self) -> int:
        if self._invalidated:
            raise ConsumedObjectError("This is an invalidated object")
        return self._num_summands

    def __iadd__(self,
                 other : Symbol[Tensorand,T]):
        """
        x += y
        for x+=x
        the id's are the same, so that is a special case
        """
        #pylint:disable=protected-access
        if self._invalidated or other._invalidated:
            raise ConsumedObjectError("This is an invalidated object")
        if id(other) == id(self):
            self._underlying = map(lambda z: (2*z[0],z[1]), self._underlying)
            return
        self._underlying = chain(self._underlying, other._underlying)
        self._num_summands += other._num_summands
        other._invalidated = True

    def __isub__(self,
                 other : Symbol[Tensorand,T]):
        if self._invalidated or other._invalidated:
            raise ConsumedObjectError("This is an invalidated object")
        if id(self) == id(other):
            self._underlying = []
            self._num_summands = 0
            return
        other.negate()
        self._underlying = chain(self._underlying, other._underlying)
        self._num_summands += other._num_summands
        other._invalidated = True

    def negate(self):
        """
        inplace negation
        """
        if self._invalidated:
            raise ConsumedObjectError("This is an invalidated object")
        self._underlying = map(lambda z: (-z[0],z[1]), self._underlying)

    def __imul__(self, scalar : float):
        if self._invalidated:
            raise ConsumedObjectError("This is an invalidated object")
        self._underlying = map(lambda z: (scalar*z[0],z[1]), self._underlying)

    def clone(self, letter_copy: Callable[[Tensorand],Tensorand]) -> Symbol[Tensorand,T]:
        """
        clone
        """
        if self._invalidated:
            raise ConsumedObjectError("This is an invalidated object")
        self._underlying = list(self._underlying)
        new_summands = [(cur_coeff,cur_summand.clone(letter_copy))
                        for (cur_coeff,cur_summand) in self._underlying]
        return Symbol(new_summands)

    #pylint:disable=too-many-arguments,too-many-locals,too-many-branches
    def ishuffle_many(self, rest: List[Symbol[Tensorand,T]],
                      materialize_length = 200,
                      do_multiprocess=True,
                      how_many_more_times_do_multiprocess = 2,
                      how_to_clone_repeats: Optional[Callable[[Tensorand],Tensorand]]=None,
                      is_top_call = True):
        """
        replace self with shuffle product of self and all the rest
        using the associativity in order to split the computation
        if the number of factors is too large
        all of the original Symbol's are consumed in the process
        TODO test
        """
        if is_top_call:
            #pylint:disable=protected-access
            if any((z._invalidated for z in rest)) or self._invalidated:
                raise ConsumedObjectError("This is an invalidated object")
            if how_to_clone_repeats is not None:
                success, rest = \
                    avoid_common_pointers(rest,
                                          lambda symbol: symbol.clone(how_to_clone_repeats))
            else:
                success, rest = avoid_common_pointers(rest, None)
            if not success:
                raise SameObjectError(
                    "The objects must have distinct ids. We don't know how to clone for you.")
            if how_to_clone_repeats is not None:
                success, rest = \
                    avoid_this_pointer(self, rest,
                                       lambda symbol: symbol.clone(how_to_clone_repeats))
            else:
                success, rest = avoid_this_pointer(self, rest, None)
            if not success:
                raise SameObjectError(
                    "The objects must have distinct ids. We don't know how to clone for you.")
        num_factors = len(rest)+1
        if num_factors<5 or not do_multiprocess:
            for other in rest:
                self.ishuffle_product(other,
                                      materialize_length=materialize_length,
                                      how_to_clone_repeats=how_to_clone_repeats)
        elif len(self) == 0:
            return
        else:
            split_half = num_factors//2
            with_self = rest[0:split_half]
            without_self_head = rest[split_half]
            without_self_tail = rest[(split_half+1):]
            if how_many_more_times_do_multiprocess > 0:
                do_multiprocess_next = True
                how_many_more_times_do_multiprocess_next = how_many_more_times_do_multiprocess-1
            else:
                do_multiprocess_next = False
                how_many_more_times_do_multiprocess_next = 0
            with Pool(2) as pool:
                both_results = [pool.apply_async(lambda main,remaining:
                                 main.ishuffle_many(remaining,
                                                    materialize_length = materialize_length,
                                                    do_multiprocess = do_multiprocess_next,
                                                    how_many_more_times_do_multiprocess =\
                                                        how_many_more_times_do_multiprocess_next,
                                                    how_to_clone_repeats = None,
                                                    is_top_call = False
                                                    ),
                                 cur_arg) for cur_arg in
                                 [(self,with_self),
                                  (without_self_head,without_self_tail)]]
                _ = [result.get() for result in both_results]
            self.ishuffle_product(without_self_head,
                                  materialize_length=materialize_length,
                                  how_to_clone_repeats=None)
        if is_top_call:
            for cur in rest:
                #pylint:disable=protected-access
                cur._invalidated = True

    def ishuffle_product(self, other: Symbol[Tensorand,T], *,
                         materialize_length: int = 200,
                         how_to_clone_repeats: Optional[Callable[[Tensorand],Tensorand]]=None):
        """
        replaces self with shuffle product of self and other
        only materializes into a list instead of lazy generator chaining
        if the length is small enough
        other is consumed in the process
        """
        #pylint:disable=protected-access
        if self._invalidated or other._invalidated:
            raise ConsumedObjectError("This is an invalidated object")
        if id(other) == id(self):
            #pylint:disable=pointless-string-statement
            """
            if the two objects were the same in memory
            and underlying was a lazy generator
            then in the double loop we would be consuming
            the same iterator in two different places
            """
            if how_to_clone_repeats is None:
                raise SameObjectError(
                    "The objects must have distinct ids. We don't know how to clone for you.")
            other = other.clone(how_to_clone_repeats)
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
        if self._num_summands<materialize_length:
            self._underlying = list(self._underlying)
        other._invalidated = True
