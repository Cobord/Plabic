"""
work arounds for two mutable references
inputs to functions thought of as having taken ownership
so anybody still with that pointer has results that can't be
relied on
e.g.
    y += x
    x is no longer a good object because of the class x,y belong too
    and the way __iadd__ was written in that class
"""

from collections import defaultdict
from typing import Callable, Dict, Iterable, List, Optional, Tuple, TypeVar

T = TypeVar("T")

def list_duplicates(seq : Iterable[T]):
    """
    which indices in the sequence were duplicates of something already seen
    """
    tally : Dict[T,List[int]] = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items()
            if len(locs)>1)

def avoid_common_pointers(seq: List[T],
                          how_to_clone: Optional[Callable[[T],T]] = None) -> Tuple[bool,List[T]]:
    """
    if any of the sequence are pointing to the same object
    do the requisite clone on the duplicates if have that function
    """
    for (_id,idces_seen) in list_duplicates((id(z) for z in seq)):
        the_unoriginals = idces_seen[1:]
        for idx in the_unoriginals:
            if how_to_clone is not None:
                seq[idx] = how_to_clone(seq[idx])
            else:
                return (False,seq)
    return (True,seq)

def avoid_this_pointer(what_to_avoid: T, seq: List[T],
                          how_to_clone: Optional[Callable[[T],T]] = None) -> Tuple[bool,List[T]]:
    """
    any in seq that has the same id as what_id, clone it if possible
    """
    the_unoriginals = (idx for idx,z in enumerate(seq) if id(z)==id(what_to_avoid))
    for idx in the_unoriginals:
        if how_to_clone is not None:
            seq[idx] = how_to_clone(seq[idx])
        else:
            return (False,seq)
    return (True,seq)

class SameObjectError(ValueError):
    """
    the two objects are exactly the same
    so lazy properties inside them won't behave correctly
    when we have a function that expects two different
    values of the same type
    """

class ConsumedObjectError(ValueError):
    """
    this object has been consumed by being the argument
    in another method, so it should be garbage collected
    already (or at least ready for garbage collection if that hasn't happened yet)
    but we don't have any ownership semantics to enforce not reusing
    so we resort to throwing an exception
    """
