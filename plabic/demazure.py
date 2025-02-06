"""
0-Hecke calculations in type A
"""

from typing import List, Optional

# does not actually check that it is a permutation
OneLinePermutation = List[int]

#pylint:disable=too-many-arguments,too-many-branches
def hopping_operator(t: int,
                     l_reversed: List[int],
                     w: OneLinePermutation,
                     n: Optional[int] = None,
                     check_inputs: bool = False,
                     where_t: Optional[int] = None) -> OneLinePermutation:
    """
    https://arxiv.org/pdf/2411.18584
    Definition 2.1
    but take reversed L instead of L
    """
    if n is None and len(w) == 0:
        raise ValueError("n must be nonzero")
    if n is None:
        n = max(w)
    if check_inputs:
        if t<=0 or t>n:
            raise ValueError(f"t must be in [1..{n}]")
        l_size = len(l_reversed)
        l_set = set(l_reversed)
        if len(l_set) != l_size:
            raise ValueError(f"L must be an ordered subset of [1..{n}]")
        for l_cur in l_reversed:
            if l_cur<=0 or l_cur>n:
                raise ValueError(f"L must be an ordered subset of [1..{n}]")
        if len(w) != n:
            raise ValueError(f"w must be a permutation of {n}")
        w_copy = w.copy()
        w_copy.sort()
        if w_copy != list(range(n)):
            raise ValueError(f"w must be a permutation of {n}")
    if where_t is None:
        idx_t = w.index(t)
    else:
        idx_t = w.index(t,where_t)
    for l_cur in l_reversed:
        if l_cur <= t:
            continue
        try:
            idx_l_cur = w.index(l_cur, idx_t)
            w[idx_t], w[idx_l_cur] = w[idx_l_cur], w[idx_t]
            return hopping_operator(t, l_reversed, w, n, False, idx_l_cur)
        except ValueError:
            pass
    return w

def up_arrow(w: OneLinePermutation,
             a: int,
             n: Optional[int] = None,
             check_input : bool = False) -> List[int]:
    """
    https://arxiv.org/pdf/2411.18584
    Definition 2.5
    but return reversed L instead
    """
    if n is None and len(w) == 0:
        raise ValueError("n must be nonzero")
    if n is None:
        n = max(w)
    if check_input:
        if a<=0 or a>n:
            raise ValueError(f"a must be in [1..{n}]")
        if len(w) != n:
            raise ValueError(f"w must be a permutation of {n}")
        w_copy = w.copy()
        w_copy.sort()
        if w_copy != list(range(n)):
            raise ValueError(f"w must be a permutation of {n}")
    idx_a = w.index(a)
    l_forward = [z for z in w[0:idx_a] if z > a]
    l_forward.reverse()
    return l_forward
