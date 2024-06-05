"""
utilities for lists being treated as items in (counter)clockwise order
"""

from itertools import chain
import itertools
from typing import List, Optional, Protocol, Tuple, TypeVar
T = TypeVar("T")

def cyclic_equal(a_list: List[T], b_list: List[T]) -> bool:
    """
    some cyclic permutation of a_list is equal to b_list
    """
    if len(a_list) != len(b_list):
        return False
    for shift in range(len(a_list)):
        up_to = len(b_list) - shift
        first_part = a_list[shift:len(a_list)] == b_list[0:up_to]
        second_part = a_list[0:shift] == b_list[up_to:len(b_list)]
        if first_part and second_part:
            return True
    return False

#pylint:disable=too-few-public-methods
class Comparable(Protocol):
    """
    the alphabet for a Lyndon word is ordered
    """
    def __leq__(self,other)->bool:...
U = TypeVar("U",bound=Comparable)

def is_lyndon(a_list: List[U]) -> bool:
    """
    is the word smaller than all of it's proper suffixes
    """
    if len(a_list)==0:
        return False
    for split in range(1,len(a_list)):
        suffix = a_list[split:]
        if not a_list<=suffix:
            return False
    return True

def generate_lyndons(ordered_alphabet : List[U], max_word_length : int):
    """
    generate the lyndon words in this ordered alphabet
    starting from shortest up to a maximum word length
    just the brute force generate and filter approach
    if this becomes a hot path then replace with Duval
    """
    num_letters = len(ordered_alphabet)
    if num_letters == 0:
        raise ValueError("The alphabet must be nonempty")
    for lyndon_length in range(1,max_word_length+1):
        for first_idx in range(num_letters):
            for rest_word in map(list,
                                 itertools.permutations(
                                     range(first_idx,num_letters),lyndon_length-1)):
                int_word = [first_idx]+rest_word
                if is_lyndon(int_word):
                    yield [ordered_alphabet[idx] for idx in int_word]

def neighboring_indices(face: List[T], this_vertex: T, that_vertex: T) -> Optional[Tuple[int,int]]:
    """
    face is a list with this_vertex and that_vertex next to each other (cyclically)
    what are their respective indices, None if that assumption not holding
    """
    if this_vertex in face and that_vertex in face:
        this_idx = face.index(this_vertex)
        that_idx = face.index(that_vertex)
        if abs(this_idx-that_idx) == 1:
            return [this_idx,that_idx]
        if this_idx==0 and that_idx==face.len()-1:
            return [this_idx,that_idx]
        if that_idx==0 and this_idx==face.len()-1:
            return [this_idx,that_idx]
        return None
    return None

def put_in_between(face : List[T], this_vertex : T, that_vertex : T, desired_name : T) -> List[T]:
    """
    face is a list with this_vertex and that_vertex next to each other (cyclically)
    put desired_name right in between
    """
    if this_vertex in face and that_vertex in face:
        this_idx = face.index(this_vertex)
        that_idx = face.index(that_vertex)
        if abs(this_idx-that_idx) == 1:
            face[this_idx] = desired_name
            _ = face.pop(that_idx)
            return face
        if this_idx==0 and that_idx==face.len()-1:
            face[this_idx] = desired_name
            _ = face.pop(that_idx)
            return face
        if that_idx==0 and this_idx==face.len()-1:
            face[that_idx] = desired_name
            _ = face.pop(this_idx)
            return face
        raise ValueError("The specified vertices weren't neighbors on face")
    raise ValueError("The specified vertices weren't on the face")

#pylint:disable=too-many-branches
def combine_cyclicly_ordered(from_self : List[T],
                             from_other : List[T],
                             glued_vertex_1 : T,
                             glued_vertex_2 : T) -> \
                                Tuple[bool, List[T], str]:
    """
    splice these two together
    """
    self_answer = neighboring_indices(from_self,glued_vertex_1,glued_vertex_2)
    if self_answer is not None:
        glued_1_idx_self, glued_2_idx_self = self_answer
        if glued_1_idx_self < glued_2_idx_self and glued_2_idx_self-glued_1_idx_self==1:
            # from_self = [...,1,2,....]
            from_self = list(chain(from_self[glued_2_idx_self+1:],from_self[0:glued_1_idx_self]))
            this_before_that_self = True
        elif glued_1_idx_self < glued_2_idx_self:
            # from_self = [1,...,2]
            from_self = list(from_self[1:len(from_self)-1])
            this_before_that_self = True
        elif glued_2_idx_self < glued_1_idx_self and glued_1_idx_self-glued_2_idx_self==1:
            # from_self = [...,2,1,....]
            from_self = list(chain(from_self[glued_1_idx_self+1:],from_self[0:glued_2_idx_self]))
            this_before_that_self = False
        elif glued_1_idx_self < glued_2_idx_self:
            # from_self = [2,...,1]
            from_self = list(from_self[1:len(from_self)-1])
            this_before_that_self = False
        else:
            raise ValueError(
                "neighboring_indices gave an answer, but it did not meet the postcriterion")
    else:
        return False, [], "said glued vertices weren't two neighbors in self's face"
    other_answer = neighboring_indices(from_other,glued_vertex_1,glued_vertex_2)
    if other_answer is not None:
        glued_1_idx_other, glued_2_idx_other = other_answer
        if glued_1_idx_other < glued_2_idx_other and \
            glued_2_idx_other-glued_1_idx_other==1:
            # from_other = [...,1,2,....]
            from_other = list(chain(from_other[glued_2_idx_other+1:],
                                    from_other[0:glued_1_idx_other]))
            this_before_that_other = True
        elif glued_1_idx_other < glued_2_idx_other:
            # from_other = [1,...,2]
            from_other = list(from_other[1:len(from_other)-1])
            this_before_that_other = True
        elif glued_2_idx_other < glued_1_idx_other and \
            glued_1_idx_other-glued_2_idx_other==1:
            # from_other = [...,2,1,....]
            from_other = list(chain(from_other[glued_1_idx_other+1:],
                                    from_other[0:glued_2_idx_other]))
            this_before_that_other = False
        elif glued_1_idx_other < glued_2_idx_other:
            # from_other = [2,...,1]
            from_other = list(from_other[1:len(from_other)-1])
            this_before_that_other = False
        else:
            raise ValueError(
                "neighboring_indices gave an answer, but it did not meet the postcriterion")
    else:
        return False, [], "said glued vertices weren't two neighbors in other's face"
    success = True
    if this_before_that_self and not this_before_that_other:
        # from_self = [....,1,2]
        # from_other = [....,2,1]
        new_face = list(chain(from_self,from_other))
    elif this_before_that_other and not this_before_that_self:
        # from_self = [....,2,1]
        # from_other = [....,1,2]
        new_face = list(chain(from_other,from_self))
    else:
        return False, [], "the order seen on the two glued elements should be opposites"
    return success, new_face, "Success"
