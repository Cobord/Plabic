"""
Manipulate sign patterns
"""
from __future__ import annotations
from enum import IntEnum
from functools import cached_property
from typing import Iterable, Protocol, TypeVar


class Reals(Protocol):
    """
    to be generic over number type
    """
    def __le__(self, other):
        ...
    @classmethod
    def zero(cls) -> Reals:
        """
        0
        """

T = TypeVar("T", bound=Reals)

class Sign(IntEnum):
    """
    +,-,0
    """
    POSITIVE = 1
    NEGATIVE = -1
    ZERO = 0

    @staticmethod
    def create(value: float) -> Sign:
        """
        sign of ordinary floats,
        rather than make you input list of Sign
        can just do ordinary float's and this will take care
        of using the enum to reduce to
        the only relevant type which is only Sign
        """
        if value<0:
            return Sign.NEGATIVE
        if value>0:
            return Sign.POSITIVE
        return Sign.ZERO

class SignPattern:
    """
    a list of signs
    """
    def __init__(self, values: Iterable[float]):
        #pylint:disable = unused-private-member
        self.__values = list(Sign.create(value) for value in values)
        _ = self.var

    @cached_property
    def var(self) -> int:
        """
        count sign changes
        ignore 0
        """
        try:
            first_positive = self.__values.index(Sign.POSITIVE)
        except ValueError:
            first_positive = len(self.__values)
        try:
            first_negative = self.__values.index(Sign.NEGATIVE)
        except ValueError:
            first_negative = len(self.__values)
        first_nonzero = min(first_positive, first_negative)
        if first_nonzero>=len(self.__values):
            raise ValueError("Must be something nonzero")
        ret_val = 0
        sign_currently = self.__values[first_nonzero]
        for next_sign in self.__values[first_nonzero+1:]:
            if next_sign == Sign.ZERO:
                continue
            if next_sign != sign_currently:
                sign_currently = next_sign
                ret_val += 1
        return ret_val

    @cached_property
    def varbar(self) -> int:
        """
        count sign changes
        0 can be 0+ or 0-
        whichever makes the count bigger
        """
        if Sign.ZERO not in self.__values:
            return self.var
        try:
            first_positive = self.__values.index(Sign.POSITIVE)
        except ValueError:
            first_positive = len(self.__values)
        try:
            first_negative = self.__values.index(Sign.NEGATIVE)
        except ValueError:
            first_negative = len(self.__values)
        first_nonzero = min(first_positive, first_negative)
        if first_nonzero>=len(self.__values):
            raise ValueError("Must be something nonzero")
        ret_val = first_nonzero
        try:
            next_positive = self.__values.index(Sign.POSITIVE,first_nonzero+1)
        except ValueError:
            next_positive = len(self.__values)
        try:
            next_negative = self.__values.index(Sign.NEGATIVE,first_nonzero+1)
        except ValueError:
            next_negative = len(self.__values)
        next_nonzero = min(next_positive,next_negative)
        if next_nonzero>=len(self.__values):
            trailing_zeros = len(self.__values) - first_nonzero - 1
            return ret_val + trailing_zeros
        prev_nonzero = first_nonzero
        while next_nonzero < len(self.__values):
            cur_sign = self.__values[prev_nonzero]
            next_sign = self.__values[next_nonzero]
            gap = next_nonzero - prev_nonzero - 1
            if next_sign != cur_sign:
                # + - -> 1
                # + 0 - -> + ? - = 1
                # + 0 0 - -> + - + - = 3
                # + 0 0 0 - -> + - + - - = 3
                ret_val += 1 + 2*(gap//2)
            if next_sign == cur_sign:
                # + + -> 0
                # + 0 + -> + - + = 2
                # + 0 0 + -> + - + + = 2
                # + 0 0 0 + -> + - + - + = 4
                ret_val += 2*((gap+1)//2)
            prev_nonzero = next_nonzero
            try:
                next_positive = self.__values.index(Sign.POSITIVE,prev_nonzero+1)
            except ValueError:
                next_positive = len(self.__values)
            try:
                next_negative = self.__values.index(Sign.NEGATIVE,prev_nonzero+1)
            except ValueError:
                next_negative = len(self.__values)
            next_nonzero = min(next_positive,next_negative)
        trailing_zeros = len(self.__values) - prev_nonzero - 1
        return ret_val + trailing_zeros
