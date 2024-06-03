"""
Shuffles
but lazily generated
"""

import itertools
from typing import List, TypeVar

def is_inverse_shuffle(p: List[int], n1: int, n2: int, check_p: bool = False) -> bool:
    """
    does this preserve the order of [0,n1) and [n1,n2)
    when p is a permutation of [0,n1+n2)
    """
    if check_p:
        if len(p) != n1+n2:
            raise ValueError("Not a permutation of [0,n1+n2)")
        p_copy = p.copy()
        p_copy.sort()
        if p_copy != list(range(n1+n2)):
            raise ValueError("Not a permutation of [0,n1+n2)")
    for (a,b) in zip(p[0:n1],p[1:n1]):
        if a>=b:
            return False
    for (a,b) in zip(p[n1:],p[n1+1:]):
        if a>=b:
            return False
    return True

def inverse_of_shuffles(n1 : int, n2 : int):
    """
    yield the inverse of shuffles
    """
    for p in itertools.permutations(range(n1+n2),n1+n2):
        if is_inverse_shuffle(p,n1,n2):
            yield p

T = TypeVar("T")

def shuffle_product(factor_1: List[T], factor_2: List[T]):
    """
    give the summands of the shuffle product
    """
    n1 = len(factor_1)
    n2 = len(factor_2)
    concatenated : List[T] = list(itertools.chain(factor_1, factor_2))
    for sigma_inverse in inverse_of_shuffles(n1,n2):
        sigma = [sigma_inverse.index(idx) for idx in range(n1+n2)]
        assert [sigma[idx] for idx in sigma_inverse] == list(range(n1+n2))
        yield [concatenated[idx] for idx in sigma]
