"""
A symbol alphabet is TODO
"""
from __future__ import annotations
from typing import Callable, Dict, FrozenSet, Generic, Iterable, List, Set, Tuple, TypeVar

Tensorand = TypeVar("Tensorand")
T = TypeVar("T")

#pylint:disable=too-few-public-methods
class SymbolWord(Generic[Tensorand,T]):
    """
    TODO
    """
    def __init__(self,
                 letters: List[Tensorand],
                 letters_variables: Callable[[Tensorand],Iterable[str]],
                 dlog : Callable[[Tensorand],Dict[str,T]]
                 ):
        if len(letters) == 0:
            raise ValueError("There should be at least one tensorand")
        self.__underlying = letters
        self.__used_letters : Set[str] = set()
        self.__no_variable_idces = []
        for (cur_idx,letter) in enumerate(letters):
            vars_in_this_letter = 0
            for var in letters_variables(letter):
                vars_in_this_letter = 0
                self.__used_letters.add(var)
            if vars_in_this_letter == 0:
                self.__no_variable_idces.append(cur_idx)
        self.__dlog = dlog

    def to_integrands(self) -> List[Dict[str,T]]:
        """
        da/a, db/b, ... each one as an element of type T
        """
        return [self.__dlog(z) for z in self.__underlying]

    def my_variables(self) -> FrozenSet[str]:
        """
        what variables are used
        """
        return frozenset(self.__used_letters)

    def am_i_zero(self) -> bool:
        """
        one of the tensorands 
        """
        return False

class Symbol:
    """
    TODO
    """
    def __init__(self, summands: List[Tuple[float,SymbolWord]]):
        #pylint:disable=unused-private-member
        self.__underlying = [(coeff,summand) for (coeff,summand)
                             in summands if not summand.am_i_zero()]

    def ishuffle_product(self, other: Symbol):
        """
        replaces self with shuffle product of self and other
        """
        raise NotImplementedError
