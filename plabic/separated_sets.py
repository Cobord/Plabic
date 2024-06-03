"""
Weakly/Strongly separated collections
"""

from abc import ABC, abstractmethod
from functools import reduce
from math import comb
from typing import List, Optional, Set, Tuple

class SeparatedCollection(ABC):
    """
    common to both strongly and weakly separated collections
    """

    def __init__(self,my_sets : List[Set[int]],*,
                 my_min : Optional[int] = None,
                 my_max : Optional[int] = None):
        self.my_sets : List[Set[int]] = []
        self.my_min = 0
        self.my_max = 0
        all_nonempty = (z.copy() for z in my_sets if len(z)>0)
        first_nonempty = next(all_nonempty)
        if first_nonempty is None:
            if my_min is None or my_max is None:
                raise ValueError("When no nonempty sets are provided, the min/max are not inferred")
            self.my_min = my_min
            self.my_max = my_max
        else:
            self.my_min = reduce(min,first_nonempty)
            self.my_max = reduce(max,first_nonempty)

    def _init_last_part(self,my_min : Optional[int],my_max : Optional[int]):
        """
        check and fix up min and max
        """
        if my_min is not None:
            if my_min > self.my_min:
                raise ValueError(f"A set had something smaller than {my_min}")
            self.my_min = min(self.my_min,my_min)
        if my_max is not None:
            if my_max < self.my_max:
                raise ValueError(f"A set had something larger than {my_max}")
            self.my_max = max(self.my_max,my_max)

    def _append_nocheck(self,next_set : Set[int]):
        """
        append another set with no checks
        """
        self.my_min = reduce(min,next_set,self.my_min)
        self.my_max = reduce(max,next_set,self.my_max)
        self.my_sets.append(next_set)

    def __len__(self):
        return len(self.my_sets)

    @abstractmethod
    def _append(self,next_set : Set[int]):
        """
        append another set making sure it is weakly/strongly separated from the rest
        for internal purposes, so it is allowed to change the range
        """

    def append(self,next_set : Set[int]):
        """
        append another set making sure it is weakly/strongly separated from the rest
        but for external purposes, so it is also not allowed to change the range
        """
        old_min = self.my_min
        old_max = self.my_max
        self._append(next_set)
        if self.my_min != old_min or self.my_max != old_max:
            bad_set = self.my_sets.pop()
            raise ValueError(f"{bad_set} was not a subset of [{old_min}, {old_max}]")

    @abstractmethod
    def maximal_by_size(self) -> bool:
        """
        how many sets can be in the collection
        """

class StronglySeparatedCollection(SeparatedCollection):
    """
    a set of sets of integers such that for any pair they are strongly separated
    """

    #pylint:disable=too-many-return-statements,too-many-branches
    @staticmethod
    def strongly_separated_disproof(i_set : Set[int], j_set : Set[int],*,
                                n_val : Optional[int] = None,
                                min_val : Optional[int] = None,
                                did_reverse = False) -> Optional[Tuple[int,int,int]]:
        """
        a disproof of the claim that i_set and j_set are strongly separated
        meaning an a<b<c such that either
        - a,c in I but not J and b in J but not I
        - a,c in J but not I and b in I but not J
        """
        if len(i_set)==0 or len(j_set)==0:
            return None
        if n_val is None:
            n_val = reduce(max,j_set,reduce(max,i_set))
        i_list = list(i_set)
        i_list.sort()
        j_list = list(j_set)
        j_list.sort()
        if i_list == list(range(n_val-len(i_set)+1,n_val+1)):
            return None
        if j_list == list(range(n_val-len(i_set)+1,n_val+1)):
            return None
        if min_val is None:
            min_val = reduce(min,j_set,reduce(min,i_set,n_val))
        if i_list == list(range(min_val,min_val+len(i_set))):
            return None
        if j_list == list(range(min_val,min_val+len(j_set))):
            return None
        for a_val in range(n_val+1):
            if a_val not in i_set or a_val in j_set:
                continue
            for b_val in range(a_val+1,n_val+1):
                if b_val not in j_set or b_val in i_set:
                    continue
                for c_val in range(b_val+1,n_val+1):
                    if c_val in i_set and c_val not in j_set:
                        return (a_val,b_val,c_val)
        if not did_reverse:
            return StronglySeparatedCollection.strongly_separated_disproof(j_set,i_set,
                                            n_val=n_val,min_val=min_val,
                                            did_reverse=True)
        return None

    def __init__(self,my_sets : List[Set[int]],*,
                 my_min : Optional[int] = None,
                 my_max : Optional[int] = None):
        super().__init__(my_sets,my_min=my_min,my_max=my_max)
        for cur_set in my_sets:
            self._append(cur_set)
        super()._init_last_part(my_min,my_max)

    def _append(self,next_set : Set[int]):
        """
        append another set making sure it is strongly separated from the rest
        internal purpose, so it does change the range
        """
        if self.maximal_by_size():
            unchanged_range = self.my_min == reduce(min,next_set,self.my_min) and \
                self.my_max == reduce(max,next_set,self.my_max)
            if unchanged_range:
                if all(z_set != next_set for z_set in self.my_sets):
                    #pylint:disable=line-too-long
                    raise ValueError(
                        f"{next_set} was not strongly separated with something because we are already at maximum {len(self)}")
                return
        ignore = False
        for prev_set in self.my_sets:
            if prev_set == next_set:
                ignore = True
                break
            next_disproof = self.strongly_separated_disproof(prev_set,next_set)
            if next_disproof is not None:
                #pylint:disable=line-too-long
                raise ValueError(
                    f"{prev_set} and {next_set} were not strongly separated because of {next_disproof}")
        if not ignore:
            super()._append_nocheck(next_set)

    def maximal_by_size(self) -> bool:
        """
        could we append another set or not
        """
        n_val = self.my_max - self.my_min + 1
        expected_num_sets = comb(n_val,0)+comb(n_val,1)+comb(n_val,2)
        return expected_num_sets == len(self)

class WeaklySeparatedCollection(SeparatedCollection):
    """
    a set of sets of integers such that for any pair they are weakly separated
    """

    #pylint:disable=too-many-return-statements,too-many-branches
    @staticmethod
    def weakly_separated_disproof(i_set : Set[int], j_set : Set[int],*,
                                n_val : Optional[int] = None,
                                min_val : Optional[int] = None,
                                did_reverse = False) -> Optional[Tuple[int,int,int,int]]:
        """
        a disproof of the claim that i_set and j_set are weakly separated
        meaning an a<b<c<d such that either
        - a,c in I but not J and b,d in J but not I
        - a,c in J but not I and b,d in I but not J
        """
        if len(i_set)==0 or len(j_set)==0:
            return None
        if n_val is None:
            n_val = reduce(max,j_set,reduce(max,i_set))
        i_list = list(i_set)
        i_list.sort()
        j_list = list(j_set)
        j_list.sort()
        if i_list == list(range(i_list[0],i_list[-1]+1)):
            return None
        if j_list == list(range(j_list[0],j_list[-1]+1)):
            return None
        if min_val is None:
            min_val = reduce(min,j_set,reduce(min,i_set,n_val))
        comp_i_set = set(range(min_val,n_val+1))
        comp_i_set.difference_update(i_list)
        comp_i_list = list(comp_i_set)
        comp_i_list.sort()
        if comp_i_list == list(range(comp_i_list[0],comp_i_list[-1]+1)):
            return None
        comp_j_set = set(range(min_val,n_val+1))
        comp_j_set.difference_update(j_list)
        comp_j_list = list(comp_j_set)
        comp_j_list.sort()
        if comp_j_list == list(range(comp_j_list[0],comp_j_list[-1]+1)):
            return None
        for a_val in range(n_val+1):
            if a_val not in i_set or a_val in j_set:
                continue
            for b_val in range(a_val+1,n_val+1):
                if b_val not in j_set or b_val in i_set:
                    continue
                for c_val in range(b_val+1,n_val+1):
                    if c_val not in i_set or c_val in j_set:
                        continue
                    for d_val in range(c_val+1,n_val+1):
                        if d_val in j_set and d_val not in i_set:
                            return (a_val,b_val,c_val,d_val)
        if not did_reverse:
            return WeaklySeparatedCollection.weakly_separated_disproof(j_set,i_set,
                                            n_val=n_val,min_val=min_val,
                                            did_reverse=True)
        return None

    #pylint:disable=too-many-arguments
    def __init__(self,my_sets : List[Set[int]],*,
                 my_min : Optional[int] = None,
                 my_max : Optional[int] = None,
                 uniform_k : bool,
                 all_sets_k : Optional[int] = None):
        super().__init__(my_sets,my_min=my_min,my_max=my_max)
        self.all_sets_k = all_sets_k
        if not uniform_k:
            all_sets_k = None
        if len(my_sets)==0 and uniform_k and all_sets_k is None:
            #pylint:disable=line-too-long
            raise ValueError(
                "When there are no sets provided and all the lengths are the same, the common set size must be provided")
        for cur_set in my_sets:
            if uniform_k and self.all_sets_k is None:
                self.all_sets_k = len(cur_set)
            self._append(cur_set)
        super()._init_last_part(my_min,my_max)

    def _append(self,next_set : Set[int]):
        """
        append another set making sure it is weakly separated from the rest
        internal purpose, so it does change the range
        """
        if self.all_sets_k is not None and\
            self.all_sets_k != len(next_set):
            raise ValueError(
                f"All the sets must have {self.all_sets_k} elements")
        if self.maximal_by_size():
            unchanged_range = self.my_min == reduce(min,next_set,self.my_min) and \
                self.my_max == reduce(max,next_set,self.my_max)
            old_set_sizes = {len(z_set) for z_set in self.my_sets}
            if unchanged_range and len(next_set) in old_set_sizes:
                if all(z_set != next_set for z_set in self.my_sets):
                    #pylint:disable=line-too-long
                    raise ValueError(
                        f"{next_set} was not weakly separated with something because we are already at maximum {len(self)}")
                return
        ignore = False
        for prev_set in self.my_sets:
            if prev_set == next_set:
                ignore = True
                break
            next_disproof = self.weakly_separated_disproof(prev_set,next_set)
            if next_disproof is not None:
                #pylint:disable=line-too-long
                raise ValueError(
                    f"{prev_set} and {next_set} were not weakly separated because of {next_disproof}")
        if not ignore:
            super()._append_nocheck(next_set)

    def maximal_by_size(self) -> bool:
        """
        could we append another set or not
        how many sets are maximal depends on whether or not
        all the sets in the collection must have some uniform size
        k_val or they can be arbitrary subsets in the range
        """
        n_val = self.my_max - self.my_min + 1
        k_val = self.all_sets_k
        if k_val is not None:
            expected_num_sets = k_val*(n_val-k_val)+1
        else:
            expected_num_sets = comb(n_val,0)+comb(n_val,1)+comb(n_val,2) + comb(n_val,3)
        return expected_num_sets == len(self)
