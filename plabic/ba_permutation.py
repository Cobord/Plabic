"""
Bounded Affine Permutations
"""

from __future__ import annotations
from copy import deepcopy
import itertools
from math import gcd
from typing import Dict, Iterator, List, Optional, Tuple, cast, Set

#pylint:disable=invalid-name
def lcm(a, b):
    """
    math only has gcd available
    """
    return abs(a*b) // gcd(a, b)

class ShiftEquivariantZBijection:
    """
    Common to all bijections of Z
    satisfying f(k+n) = f(k)+n for all k
    for some n
    """
    my_n : int
    my_f : Dict[int,int]
    is_lift_from_sn : bool

    # pylint:disable=too-many-locals
    def __init__(self,*,
                 coxeter_word : Optional[List[int]]=None,
                 some_vals : Optional[Dict[int,int]]=None,
                 juggling_pattern : Optional[List[int]]=None,
                 my_sn : Optional[List[int]]=None,
                 decorations : Optional[List[Optional[bool]]]=None,
                 to_bounded_affine : Optional[bool] = None,
                 n_val : Optional[int] = None,
                 check_multiple_constructions : bool = True):
        construction_1 = my_sn is not None
        construction_2 = some_vals is not None
        construction_3 = juggling_pattern is not None
        construction_4 = coxeter_word is not None and n_val is not None
        how_many_of_123 = \
            sum(1 if z_var else 0 for z_var in [construction_1,construction_3,construction_4])
        if how_many_of_123>1:
            raise ValueError("Only one option should be used otherwise we would have to"+\
                             "test if they gave the same results")
        used_construction_4 = False
        if construction_1:
            my_sn = cast(List[int],my_sn)
            if decorations is None:
                decorations = [None]*len(my_sn)
            if to_bounded_affine is None:
                to_bounded_affine = False
            self.__lift_from_sn(my_sn,
                                decorations,
                                to_bounded_affine)
        elif construction_2:
            some_vals = cast(Dict[int,int],some_vals)
            self.__from_some_values(some_vals,n_val)
        elif construction_3:
            juggling_pattern = cast(List[int],juggling_pattern)
            self.__from_juggling_pattern(juggling_pattern)
        elif construction_4:
            coxeter_word = cast(List[int],coxeter_word)
            n_val = cast(int,n_val)
            self.__from_coxeter_word(n_val,coxeter_word)
            used_construction_4 = True
        else:
            raise ValueError("Must provide some way of construction")
        if check_multiple_constructions and not used_construction_4 and construction_4:
            from_coxeter = ShiftEquivariantZBijection(coxeter_word=coxeter_word,n_val=n_val)
            if self != from_coxeter:
                raise ValueError("The Coxeter word was inconsistent with other methods"+\
                                 "of specifying this ShiftEquivariantZBijection")
        self.is_lift_from_sn = all(1<=self[idx]<=self.my_n for idx in range(1,self.my_n+1))

    # pylint:disable=too-many-locals
    def __from_some_values(self,some_vals : Dict[int,int], n_val : Optional[int] = None):
        if n_val is None:
            n_val = len(some_vals)
        if len(some_vals)<n_val:
            raise ValueError(" ".join([
                "Not enough values of f provided",
                "to make a ShiftEquivariantZBijection",
                f"with shift {n_val}"]))
        self.my_n = n_val
        self.my_f : Dict[int,int] = {}

        periodicity_error = ValueError(" ".join([
                    "This did not meet the condition",
                    "of having shift equivariance",
                    f"with shift {self.my_n}"]))

        for (key1,value1),(key2,value2) in itertools.combinations(some_vals.items(),2):
            _,key_diff_rem = divmod(key1 - key2,n_val)
            if key_diff_rem == 0 and value1-value2 != key1-key2:
                raise ValueError(" ".join([
                    f"Specified values for {key1} and {key2}",
                    f"but they did not meet the {n_val} shift equivariance"]))
        for (key,value) in some_vals.items():
            (quot,rem) = divmod(key,self.my_n)
            if not self.__try_to_add_value(rem,value-quot*self.my_n):
                raise periodicity_error
            if rem<=0:
                if not self.__try_to_add_value(rem+self.my_n,self.my_f[rem]+self.my_n):
                    raise periodicity_error
            if not self.__try_to_add_value(key,value):
                raise periodicity_error
        for idx in range(1,self.my_n+1):
            quot_f_idx = self[idx] % self.my_n
            for jdx in range(idx+1,self.my_n+1):
                quot_f_jdx = self[jdx] % self.my_n
                if quot_f_jdx == quot_f_idx:
                    raise ValueError("This does not define a self-bijection of the integers")

    def __lift_from_sn(self,my_sn : List[int],
                     decorations : List[Optional[bool]],
                     to_bounded_affine : bool):
        """
        Ends up either an affine permutation or a bounded affine permutation
        - If an affine permutation, we know where we want to send 1-n from my_sn
        - If a (k,n) bounded affine permutation, we need the decorated permutation
            where k is determined by the decorations at fixed points
            The decoration says whether or not it is a coloop in terms of Plabic Graph
        """
        my_n = len(my_sn)
        if sorted(my_sn) != list(range(1,my_n+1)):
            raise ValueError(f"The first argument had to be in S_{my_n}")
        self.my_n = my_n
        if to_bounded_affine and len(decorations)<my_n:
            decorations += [None]*my_n
        input_dict : Dict[int,int] = {}
        for (idx,(value,decoration)) in enumerate(zip(my_sn,decorations)):
            if not to_bounded_affine:
                input_dict[idx+1] = value
            elif idx+1==value:
                if decoration is None:
                    raise ValueError("Fixed points must have decorations")
                input_dict[idx+1] = value + (my_n if decoration else 0)
            elif idx+1<value:
                input_dict[idx+1] = value
            else:
                input_dict[idx+1] = value + my_n
        self.my_f = input_dict

    def __from_coxeter_word(self,my_n : int,cox_word : List[int]):
        """
        from a word in the coxeter generators s_0 through s_(n-1)
        ends up being an Affine Permutation
        """
        self.my_n = my_n
        input_dict = {z:z for z in range(1,my_n+1)}
        for letter in cox_word:
            letter = letter % my_n
            if 1<=letter<my_n:
                (input_dict[letter+1], input_dict[letter]) = \
                    (input_dict[letter], input_dict[letter+1])
            else:
                (input_dict[1], input_dict[my_n]) = \
                    (input_dict[my_n]-my_n, input_dict[1]+my_n)
        self.my_f = input_dict

    def __from_juggling_pattern(self,juggling_pattern : List[int]):
        """
        from a juggling pattern
        """
        my_n = len(juggling_pattern)
        average_of_juggle = sum(juggling_pattern)
        if average_of_juggle % my_n != 0:
            raise ValueError("Not a valid juggling pattern")
        average_of_juggle //= my_n
        ehrenborg_readdy : Dict[int,int] = {}
        for idx in range(1,my_n+1):
            ehrenborg_readdy[idx] = idx + juggling_pattern[idx-1] - average_of_juggle
        try:
            self.__from_some_values(some_vals=ehrenborg_readdy,n_val=my_n)
        except ValueError as exc:
            raise ValueError("Not a valid juggling pattern to define"+\
                             "a shift equivariant bijection of Z") from exc

    def __try_to_add_value(self,key : int,value : int) -> bool:
        """
        try to add this key value pair to my_f
        if there was something different already there return False
        so the caller can raise the proper Exception
        """
        old_val = self.my_f.get(key,None)
        if old_val is None:
            self.my_f[key] = value
        elif old_val != value:
            return False
        return True

    def __getitem__(self,key : int) -> int:
        if (key_value := self.my_f.get(key,None)) is not None:
            return key_value
        (quot,rem) = divmod(key,self.my_n)
        rem_value = self.my_f.get(rem,None)
        if rem_value is None:
            rem_value = self.my_f.get(rem+self.my_n,None)
            if rem_value is not None:
                rem_value -= self.my_n
            else:
                raise ValueError(f"Didn't have the value for {rem} or {rem+self.my_n}")
            self.__try_to_add_value(rem,rem_value)
        key_value = rem_value+quot*self.my_n
        return key_value

    def __eq__(self,other) -> bool:
        """
        are they the same shift-equivariant bijection
        with the same shift equivariance
        """
        if not isinstance(other,ShiftEquivariantZBijection):
            return False
        if self.my_n != other.my_n:
            return False
        self_window = [self[z] for z in range(1,self.my_n+1)]
        other_window = [other[z] for z in range(1,other.my_n+1)]
        return self_window == other_window

    def __hash__(self) -> int:
        self_window = [self[z] for z in range(1,self.my_n+1)]
        return hash(" ".join(str(z_var) for z_var in self_window))

    def __mul__(self,other : ShiftEquivariantZBijection) -> ShiftEquivariantZBijection:
        """
        Do other first then self
        to get the action of self*other
        Then use those points to make an ShiftEquivariantZBijection
        """
        final_n = lcm(self.my_n,other.my_n)
        return ShiftEquivariantZBijection(
            some_vals={z : self[other[z]] for z in range(1,final_n+1)},
            n_val=final_n)

    def quot_to_sn(self) -> List[int]:
        """
        Forget down to the residues modulo self.my_n
        """
        def one_through_n_mod(some_num : int) -> int:
            if self.is_lift_from_sn:
                return some_num
            rem = some_num % self.my_n
            if rem == 0:
                rem += self.my_n
            return rem
        return [one_through_n_mod(self[z]) for z in range(1,self.my_n+1)]

    def __str__(self) -> str:
        zeroth_line = f"{self.my_n} shift equivariant bijection of Z"
        first_line = ",".join((str(i) for i in range(1,self.my_n+1)))
        second_line = ",".join((str(self[i]) for i in range(1,self.my_n+1)))
        return f"{zeroth_line}\n{first_line}\n{second_line}"

class BoundedAffinePermutation:
    """
    A (k,n) bounded affine permutation
    """
    my_underlying : ShiftEquivariantZBijection
    my_k : int
    my_n : int

    def __init__(self,*,
                 coxeter_word : Optional[List[int]]=None,
                 some_vals : Optional[Dict[int,int]]=None,
                 juggling_pattern : Optional[List[int]]=None,
                 my_sn : Optional[List[int]]=None,
                 decorations : Optional[List[Optional[bool]]]=None,
                 n_val : Optional[int] = None,
                 check_multiple_constructions : bool = True):
        self.my_underlying = ShiftEquivariantZBijection(coxeter_word=coxeter_word,
                 some_vals=some_vals,
                 juggling_pattern=juggling_pattern,
                 my_sn=my_sn,
                 decorations=decorations,
                 to_bounded_affine=True,
                 n_val=n_val,check_multiple_constructions=check_multiple_constructions)
        self.my_n = self.my_underlying.my_n

        my_invariant = sum(self[idx]-idx for idx in range(1,self.my_n+1))
        (my_k,rem) = divmod(my_invariant,self.my_n)
        # this assertion is covered by the previous checks for shift equivariance and self-bijection
        assert rem==0, \
            f"The sum of f(i)-i as i goes from 1 through n should be a multiple of {self.my_n}"
        self.my_k = my_k
        boundedness = all(idx<=self[idx]<=idx+self.my_n for idx in range(1,self.my_n+1))
        if not boundedness:
            raise ValueError(
                " ".join(["Failed to satisfy the boundedness condition for a",
                f"({abs(self.my_k)},{self.my_n}) bounded affine permutation"]))

    def __getitem__(self, idx: int) -> int:
        return self.my_underlying[idx]

    def quot_to_sn(self) -> List[int]:
        """
        Forget down to the residues modulo self.my_n
        """
        return self.my_underlying.quot_to_sn()

    def __str__(self) -> str:
        zeroth_line = f"({self.my_k},{self.my_n}) bounded affine permutation"
        first_line = ",".join((str(i) for i in range(1,self.my_n+1)))
        second_line = ",".join((str(self[i]) for i in range(1,self.my_n+1)))
        return f"{zeroth_line}\n{first_line}\n{second_line}"

    def __eq__(self,other) -> bool:
        """
        are they the same bounded affine permutation
        with the same k and n
        """
        if not isinstance(other,BoundedAffinePermutation):
            return False
        if self.my_n != other.my_n:
            return False
        if self.my_k != other.my_k:
            return False
        return self.my_underlying == other.my_underlying

    def __hash__(self) -> int:
        self_window = [self[z] for z in range(1,self.my_n+1)]
        return hash(" ".join(str(z_var) for z_var in self_window))

    @staticmethod
    def all_bounded_affine_perms(my_n : int) -> Iterator[BoundedAffinePermutation]:
        """
        Iterate through all bounded affine permuations with shift equivariance
        given by my_n
        """
        for perm in itertools.permutations(range(1,my_n+1)):
            fixed_pts = [idx+1 for idx in range(my_n) if perm[idx]==idx+1]
            num_fixed_pts = len(fixed_pts)
            real_decoration : List[Optional[bool]] = [None]*my_n
            for my_decs in itertools.product([True,False],repeat=num_fixed_pts):
                for (cur_fixed_pt,cur_dec) in zip(fixed_pts,my_decs):
                    real_decoration[cur_fixed_pt-1] = cur_dec
                yield BoundedAffinePermutation(my_sn=list(perm),
                                               decorations=real_decoration,n_val=my_n)
                real_decoration = [None]*my_n

class AffinePermutation:
    """
    An affine permutation with some periodicity n
    """
    is_bounded : bool
    my_coxeter_word : Optional[List[int]]
    my_n : int
    is_lift_from_sn : bool

    def __init__(self,*,
                 coxeter_word : Optional[List[int]]=None,
                 some_vals : Optional[Dict[int,int]]=None,
                 juggling_pattern : Optional[List[int]]=None,
                 my_sn : Optional[List[int]]=None,
                 n_val : Optional[int] = None,
                 check_multiple_constructions : bool = True):
        self.my_underlying = ShiftEquivariantZBijection(coxeter_word=coxeter_word,
                 some_vals=some_vals,
                 juggling_pattern=juggling_pattern,
                 my_sn=my_sn,
                 decorations=None,
                 to_bounded_affine=False,
                 n_val=n_val,check_multiple_constructions=check_multiple_constructions)
        self.my_n = self.my_underlying.my_n
        if self.my_n <= 1:
            raise ValueError(f"n must be at least 2 not {self.my_n}")
        self.is_lift_from_sn = self.my_underlying.is_lift_from_sn
        if coxeter_word is not None:
            self.my_coxeter_word = coxeter_word.copy()
        else:
            self.my_coxeter_word = None

        my_invariant = sum(self[idx]-idx for idx in range(1,self.my_n+1))
        assert my_invariant==0, \
            f"The sum of f(i)-i as i goes from 1 through {self.my_n} should be 0"

    def inv(self) -> AffinePermutation:
        """
        the inverse in the group
        """
        if self.my_coxeter_word is not None:
            flipped_word = self.my_coxeter_word.copy()
            flipped_word.reverse()
            return AffinePermutation(
                 coxeter_word = flipped_word,
                 n_val = self.my_n)
        some_vals : Dict[int,int] = {}
        for idx in range(1,self.my_n+1):
            f_idx = self[idx]
            some_vals[f_idx] = idx
        return AffinePermutation(
                 some_vals = some_vals,
                 n_val = self.my_n)

    def __mul__(self,other : AffinePermutation) -> AffinePermutation:
        """
        Do other first then self
        to get the action of self*other
        Then use those points to make an AffinePermutaion
        """
        final_n = lcm(self.my_n,other.my_n)
        if self.my_coxeter_word is not None and other.my_coxeter_word is not None:
            my_new_coxeter_word = self.my_coxeter_word.copy() + other.my_coxeter_word.copy()
        else:
            my_new_coxeter_word = None
        return AffinePermutation(some_vals={z : self[other[z]] for z in range(1,final_n+1)},
                                n_val=final_n,coxeter_word=my_new_coxeter_word,
                                check_multiple_constructions=False)

    def side_by_side(self,other : AffinePermutation) -> AffinePermutation:
        """
        given two finite permutations make
        the finite permutation given by taking
        - the 1-n window of self
        - the 1-m window of other
        and concatenating them as 1-(n+m) window
        where the first n are from self unchanged
        and the remaining are from other shifted up by n
        """
        if not self.is_lift_from_sn or not other.is_lift_from_sn:
            raise ValueError("Both must be finite permutations")
        final_n = self.my_n + other.my_n
        if self.my_coxeter_word is not None and \
            other.my_coxeter_word is not None and \
            0 not in self.my_coxeter_word and \
            0 not in other.my_coxeter_word:
            my_new_coxeter_word = other.my_coxeter_word.copy()
            my_new_coxeter_word = [z_var+self.my_n for z_var in my_new_coxeter_word]
            my_new_coxeter_word = self.my_coxeter_word.copy() + my_new_coxeter_word
        else:
            my_new_coxeter_word = None
        some_vals = {}
        for idx in range(1,self.my_n+1):
            some_vals[idx] = self[idx]
        for idx in range(self.my_n+1,final_n+1):
            some_vals[idx] = other[idx-self.my_n]+self.my_n
        return AffinePermutation(some_vals=some_vals,n_val=final_n,
                                 coxeter_word=my_new_coxeter_word,
                                 check_multiple_constructions=False)

    def is_fully_commutative(self) -> bool:
        """
        Can every reduced expression in the coxeter generators for this affine permutation
        be transformed to any other using only commutation moves
        """
        for (idx,jdx,kdx) in itertools.combinations(range(1,self.my_n*3+1),3):
            if self[idx]>self[jdx]>self[kdx]:
                return False
        return True

    def __eq__(self,other):
        if not isinstance(other,AffinePermutation):
            return False
        windows_equal = self.my_underlying == other.my_underlying
        if windows_equal and self.my_coxeter_word is None and \
            other.my_coxeter_word is not None:
            self.my_coxeter_word = other.my_coxeter_word.copy()
        if windows_equal and other.my_coxeter_word is None and \
            self.my_coxeter_word is not None:
            other.my_coxeter_word = self.my_coxeter_word.copy()
        if windows_equal and self.my_coxeter_word is not None and \
            other.my_coxeter_word is not None:
            self_has_smaller = len(self.my_coxeter_word)<len(other.my_coxeter_word)
            other_has_smaller = len(other.my_coxeter_word)<len(self.my_coxeter_word)
            if self_has_smaller:
                other.my_coxeter_word = self.my_coxeter_word.copy()
            if other_has_smaller:
                self.my_coxeter_word = other.my_coxeter_word.copy()
        return windows_equal

    def __hash__(self) -> int:
        self_window = [self[z] for z in range(1,self.my_n+1)]
        return hash(" ".join(str(z_var) for z_var in self_window))

    def __str__(self) -> str:
        affineness = "" if self.is_lift_from_sn else "Affine "
        periodicity = "" if self.is_lift_from_sn else "periodicity "
        zeroth_line = f"{affineness}Permutation of {periodicity}{self.my_n}"
        first_line = ",".join((str(i) for i in range(1,self.my_n+1)))
        second_line = ",".join((str(self[i]) for i in range(1,self.my_n+1)))
        return f"{zeroth_line}\n{first_line}\n{second_line}"

    def __getitem__(self, idx: int) -> int:
        return self.my_underlying[idx]

    def bruhat_leq(self,other:AffinePermutation) -> bool:
        """
        self <= other in the Bruhat partial order
        """
        raise NotImplementedError

    def quot_to_sn(self) -> List[int]:
        """
        Forget down to the residues modulo self.my_n
        """
        return self.my_underlying.quot_to_sn()

    def __len__(self) -> int:
        inversion_count = 0
        for idx in range(1,self.my_n+1):
            for jdx in range(idx+1,self.my_n+1):
                if self[idx]>self[jdx]:
                    inversion_count += 1
            cur_batch_start = self.my_n+1
            while True:
                if all(self[jdx]>self[idx] for
                       jdx in range(cur_batch_start,cur_batch_start+self.my_n)):
                    break
                inversion_count += sum(self[jdx]<self[idx] for
                       jdx in range(cur_batch_start,cur_batch_start+self.my_n))
                cur_batch_start += self.my_n
        return inversion_count

    def right_descents(self) -> List[int]:
        """
        The right descents of this affine permutation
        these are the i such that len(self*s_i)<len(self)
        """
        return [potential for potential
                in range(0,self.my_n)
                if self[potential]>self[potential+1]]

    def is_k_grassmannian(self) -> Tuple[bool,int]:
        """
        is this k-Grassmannian for some k and if so which k
        """
        my_descents = self.right_descents()
        if len(my_descents)==1:
            return (True,my_descents[0])
        return (False,0)

    @staticmethod
    def all_k_parabolic_perms(my_n : int, my_k : int, max_length : Optional[int])\
        -> Iterator[Tuple[int,AffinePermutation]]:
        """
        all elements of S_k times S_{n-k} subset S_n subset hat S_n
        that have length below a given maximum (or the maximum possible)
        they do not come out sorted by length, but the lengths are also yielded
        """
        n_minus_k = my_n - my_k
        if max_length is None:
            max_length = my_k*(my_k-1)//2 + n_minus_k*(n_minus_k-1)//2
        part_1_possibilities = AffinePermutation.all_finite_perms(my_k,max_length)
        part_2_possibilities = AffinePermutation.all_finite_perms(n_minus_k,max_length)
        for (part_1_len,part_1_perm),(part_2_len,part_2_perm) in\
            itertools.product(part_1_possibilities,part_2_possibilities):
            combined_length = part_1_len+part_2_len
            if combined_length<=max_length:
                yield combined_length,part_1_perm.side_by_side(part_2_perm)

    @staticmethod
    def all_finite_perms(my_n : int, max_length : Optional[int] = None)\
        -> Iterator[Tuple[int,AffinePermutation]]:
        """
        Iterate through all finite permuations of my_n
        sorted by length
        """
        if max_length is None:
            max_length = my_n*(my_n-1)//2
        cox_gen = [AffinePermutation(coxeter_word=[idx],n_val=my_n) for idx in range(1,my_n)]
        prev_options = [AffinePermutation(coxeter_word=[],n_val=my_n)]
        length_emitting = 0
        if length_emitting<=max_length:
            yield length_emitting,prev_options[0]
        while True:
            length_emitting += 1
            if length_emitting>max_length:
                break
            cur_round_emitted : Set[AffinePermutation] = set([])
            for cur_prev in prev_options:
                dont_use = cur_prev.right_descents()
                potential_adds = [idx for idx in range(1,my_n) if idx not in dont_use]
                for new_gen in potential_adds:
                    might_work = cur_prev*cox_gen[new_gen-1]
                    if might_work not in cur_round_emitted:
                        cur_round_emitted.add(might_work)
                        yield length_emitting,deepcopy(might_work)
            prev_options = list(cur_round_emitted)

    @staticmethod
    def all_affine_perms(my_n : int, max_length : int)\
        -> Iterator[Tuple[int,AffinePermutation]]:
        """
        Iterate through all affine permuations with shift equivariance
        given by my_n by length
        """
        cox_gen = [AffinePermutation(coxeter_word=[idx],n_val=my_n) for idx in range(my_n)]
        prev_options = [AffinePermutation(coxeter_word=[],n_val=my_n)]
        length_emitting = 0
        if length_emitting<=max_length:
            yield length_emitting,prev_options[0]
        while True:
            length_emitting += 1
            if length_emitting>max_length:
                break
            cur_round_emitted : Set[AffinePermutation] = set([])
            for cur_prev in prev_options:
                dont_use = cur_prev.right_descents()
                potential_adds = [idx for idx in range(my_n) if idx not in dont_use]
                for new_gen in potential_adds:
                    might_work = cur_prev*cox_gen[new_gen]
                    if might_work not in cur_round_emitted:
                        cur_round_emitted.add(might_work)
                        yield length_emitting,deepcopy(might_work)
            prev_options = list(cur_round_emitted)

#pylint:disable=too-many-instance-attributes
class BruhatInterval:
    """
    An interval in the Bruhat order
    Representing the set that is
    greater than or equal to v
    less than or equal to w
    the interval is empty if and only if not v<=w
    that failure of inequality can be due to incomparibility
    rather than the reverse strict inequality
    """
    my_n : int
    v : AffinePermutation
    w : AffinePermutation
    v_len : int
    w_len : int
    len_diff : int
    v_finite : bool
    w_finite : bool
    manifestly_empty : bool

    def __init__(self,my_v : AffinePermutation, my_w : AffinePermutation):
        self.my_n = my_v.my_n
        assert my_v.my_n == my_w.my_n
        self.v = my_v
        self.w = my_w
        self.v_len = len(my_v)
        self.w_len = len(my_w)
        self.manifestly_empty = self.v_len>self.w_len
        self.len_diff = self.w_len - self.v_len
        self.v_finite = my_v.is_lift_from_sn
        self.w_finite = my_w.is_lift_from_sn

    def more_in_equivalence_class(self,k_val : int) -> Iterator[BruhatInterval]:
        """
        from interval [v,w]
        consider all intervals [vz,wz] where z is drawn from k parabolics
        also assure that the Bruhat lengths of both sides shifted by the same amount
        """
        for _, parabolic in AffinePermutation.all_k_parabolic_perms(self.my_n,k_val,None):
            temp_v = deepcopy(self.v)
            temp_w = deepcopy(self.w)
            temp_v *= parabolic
            temp_w *= parabolic
            potential_interval = BruhatInterval(temp_v,temp_w)
            if self.len_diff == potential_interval.len_diff:
                yield potential_interval

    def containment(self,other : BruhatInterval) -> bool:
        """
        is self contained within other
        """
        if self.manifestly_empty:
            return True
        actually_nonempty = self.v.bruhat_leq(self.w)
        if not actually_nonempty:
            self.manifestly_empty = True
            return True
        if other.manifestly_empty:
            return False
        actually_nonempty = other.v.bruhat_leq(other.w)
        if not actually_nonempty:
            other.manifestly_empty = True
            return False
        return other.v.bruhat_leq(self.v) and self.w.bruhat_leq(other.w)

    @staticmethod
    def all_Qkn(my_k : int,
                my_n : int,
                max_length : Optional[int] = None) -> Iterator[BruhatInterval]:
        """
        Q_kn is defined as (v,w) in S_N times S_N^k
        satisfying v<=w and that S_N^k means that it is k-Grassmannian
        """
        if max_length is None:
            max_length = my_n*(my_n-1)//2
        for w_potential_len,w_potential in AffinePermutation.all_finite_perms(my_n,max_length):
            is_k_grassmannian, which_k = w_potential.is_k_grassmannian()
            if not is_k_grassmannian or which_k != my_k:
                continue
            for _,v_potential in AffinePermutation.all_finite_perms(my_n,w_potential_len):
                if v_potential.bruhat_leq(w_potential):
                    yield BruhatInterval(v_potential,w_potential)
