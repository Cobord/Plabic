"""
A R^{k(n-k)}_{gt 0} chart in TNN Grassmannian
"""
from __future__ import annotations
from copy import deepcopy
from functools import reduce
import itertools
from math import comb
from typing import Callable, Dict, FrozenSet, Iterable, List, Optional,Set, Union
from sympy import Expr,Symbol,Integer
import numpy as np
from .sym_utils import nn_with_my_symbols,determinant

# up to caller to put this into some data structure that does convex hull
# manipulations, here it just says to take the convex hull of the provided points
ConvexHull = List[np.ndarray]

def create_convex_hull(points: Iterable[np.ndarray]) -> ConvexHull:
    """
    Creates a convex hull
    TODO a representation better for manipulation than just saying
    it is the convex hull of these points
    """
    return list(points)

#pylint:disable=too-many-instance-attributes
class MinorSignConstraints:
    """
    for each minor described by a set of int
    is it either
    - manifestly positive
    - manifestly negative
    - manifestly zero
    - manifestly nonpositive
    - manifestly nonnegative
    - manifestly nonzero
    if a possible minor is not in any of these, then it is unconstrained
    """
    #pylint:disable=too-many-arguments
    def __init__(self,
                 which_positive : List[Set[int]],
                 which_negative : List[Set[int]],
                 which_zero : List[Set[int]],
                 which_nonpositive : List[Set[int]],
                 which_nonnegative : List[Set[int]],
                 which_nonzero : List[Set[int]]):
        self.set_size = None
        self.n_at_least = 0
        for cur_list in [which_positive,which_negative,
                         which_zero,which_nonpositive,
                         which_nonnegative,which_nonzero]:
            for cur_set in cur_list:
                #pylint:disable=nested-min-max
                self.n_at_least = max(self.n_at_least,max(cur_set))
                if self.set_size is None:
                    self.set_size = len(cur_set)
                elif len(cur_set) != self.set_size:
                    raise ValueError("All sets must be the same size")
                if any(z<0 for z in cur_set):
                    raise ValueError("All sets must be of natural numbers")
        self.which_positive = which_positive
        self.which_negative = which_negative
        self.which_zero = which_zero
        self.which_nonnegative = which_nonnegative
        self.which_nonpositive = which_nonpositive
        self.which_nonzero = which_nonzero

    @staticmethod
    def all_positive(my_n : int,my_k : int) -> MinorSignConstraints:
        """
        the manifestly positive constraint on all the possibilities
        for all selection of k columns from n possible
        """
        which_positive = [set(z) for z in itertools.combinations(range(my_n),my_k)]
        return MinorSignConstraints(which_positive,[],[],[],[],[])

    @staticmethod
    def all_nonnegative(my_n : int,my_k : int) -> MinorSignConstraints:
        """
        the manifestly nonnegative constraint on all the possibilities
        for all selection of k columns from n possible
        """
        which_nonnegative = [set(z) for z in itertools.combinations(range(my_n),my_k)]
        return MinorSignConstraints([],[],[],[],which_nonnegative,[])

    def is_positroidy(self, num_cols : int) -> bool:
        """
        is it the kind of sign constraint for a positroid cell in Gr(self.set_size,num_cols)
        """
        if self.n_at_least>=num_cols:
            raise ValueError("There were sign constraints with columns that do not exist")
        if self.set_size is None:
            return False
        total_pluckers = comb(num_cols,self.set_size)
        return total_pluckers == len(self.which_positive) + len(self.which_zero)

    def which_schubert_cell(self, num_cols : int,
                            sorter : Optional[Callable[[int],int]] = None)\
                                -> Set[int]:
        """
        give the corresponding schubert cell
        possibly after using a different lexicographic order
        if you do with all n! possibilities for sorter, then this is the
        Gelfand-Serganova theorem
        """
        if self.n_at_least>=num_cols:
            raise ValueError("There were sign constraints with columns that do not exist")
        if self.set_size is None:
            raise ValueError("No sets were provided so the set size property was never set")
        total_pluckers = comb(num_cols,self.set_size)
        num_nonzeros = len(self.which_positive) + len(self.which_negative) + len(self.which_nonzero)
        if total_pluckers != len(self.which_zero) + num_nonzeros:
            raise ValueError(" ".join(["The constraints should be like a matroid specifying",
                                       "some (but not all) plucker coordinates as zero",
                                       "and the rest as nonzero"]))
        if num_nonzeros == 0:
            raise ValueError(" ".join(["The constraints should be like a matroid specifying",
                                       "some (but not all) plucker coordinates as zero",
                                       "and the rest as nonzero"]))
        if sorter is None:
            # pylint: disable=unnecessary-lambda-assignment
            sorter = lambda z : z
        def my_set_leq(set_a : Set[int], set_b : Set[int]) -> bool:
            """
            set_a sorted compared to set_b sorted as lists
            """
            list_a = list(set_a)
            list_a.sort(key = sorter)
            list_b = list(set_b)
            list_b.sort(key = sorter)
            return list_a <= list_b
        lex_min_set = reduce(lambda acc,x: x.copy() if my_set_leq(x,acc) else acc,
                             itertools.chain(self.which_positive,
                                             self.which_negative,
                                             self.which_nonzero))
        return lex_min_set

    def to_hypersimplex_positroid_polytope(self,num_cols : int) -> Optional[ConvexHull]:
        """
        if this is positroidy, give the positroid polytope in the hypersimplex
        given by the image of the moment map from this positroid in tnn grassmanian
        to the hypersimplex
        See https://arxiv.org/abs/2104.08254
        """
        if not self.is_positroidy(num_cols):
            return None
        #pylint:disable=unused-variable
        def indicator_vector(which_is : Set[int]) -> np.ndarray:
            to_return = np.zeros(num_cols)
            for cur_i in which_is:
                to_return[cur_i] = 1
            return to_return
        all_pts = (indicator_vector(z) for z in self.which_positive)
        return create_convex_hull(all_pts)

    def cyclic_shift(self,my_n : int,num_times : int) -> None:
        """
        i->i-1
        0->n-1
        in all the sets
        """
        if my_n < self.n_at_least:
            raise ValueError("There were sign constraints with columns that do not exist")
        for cur_list in [self.which_positive,self.which_negative,
                         self.which_zero,self.which_nonpositive,
                         self.which_nonnegative,self.which_nonzero]:
            for idx_in_list,cur_set in enumerate(cur_list):
                cur_list[idx_in_list] = {(z-num_times)%my_n for z in cur_set}

    #pylint:disable=too-many-return-statements
    def check(self,k_subset : Set[int], cur_expr : Expr, var_set : Set[Symbol]) -> bool:
        """
        whichever constraint k_subset is in, check that sign constraint is satisfied
        if it is in none, then there is no constraint
        """
        if k_subset in self.which_positive:
            return nn_with_my_symbols(cur_expr,True,var_set)
        if k_subset in self.which_negative:
            return nn_with_my_symbols(-cur_expr,True,var_set)
        if k_subset in self.which_zero:
            return cur_expr == Integer(0)
        if k_subset in self.which_nonnegative:
            return nn_with_my_symbols(cur_expr,False,var_set)
        if k_subset in self.which_nonpositive:
            return nn_with_my_symbols(-cur_expr,False,var_set)
        if k_subset in self.which_nonzero:
            is_positive = nn_with_my_symbols(cur_expr,True,var_set)
            is_negative = nn_with_my_symbols(-cur_expr,True,var_set)
            return is_positive or is_negative
        return True

    def to_manifestly_positroidy(self,num_cols : int) -> PositroidySignConstraints:
        """
        convert to PositroidySignConstraints
        """
        assert self.is_positroidy(num_cols)
        return PositroidySignConstraints(
            which_positive = self.which_positive,
            which_zero=self.which_zero)

class PositroidySignConstraints:
    """
    when manifestly positroidy,
    only store one of which_positive and which_zero
    the others must be empty
    and the one not being stored is all the other possibilities
    """
    def __init__(self,*,which_positive : List[Set[int]],which_zero : List[Set[int]]):
        self.set_size = None
        self.n_at_least = 0
        for cur_list in [which_positive,
                         which_zero]:
            for cur_set in cur_list:
                #pylint:disable=nested-min-max
                self.n_at_least = max(self.n_at_least,max(cur_set))
                if self.set_size is None:
                    self.set_size = len(cur_set)
                elif len(cur_set) != self.set_size:
                    raise ValueError("All sets must be the same size")
                if any(z<0 for z in cur_set):
                    raise ValueError("All sets must be of natural numbers")
        if len(which_positive)<len(which_zero):
            self.which_positive = which_positive
            self.which_zero = None
        else:
            self.which_positive = None
            self.which_zero = which_zero

    def is_positroidy(self, _num_cols : int) -> bool:
        """
        this is automatically positroidy
        """
        return True

    def which_schubert_cell(self, num_cols : int,
                            sorter : Optional[Callable[[int],int]] = None)\
                                -> Set[int]:
        """
        give the corresponding schubert cell
        possibly after using a different lexicographic order
        if you do with all n! possibilities for sorter, then this is the
        Gelfand-Serganova theorem
        """
        if self.n_at_least>=num_cols:
            raise ValueError("There were sign constraints with columns that do not exist")
        if self.set_size is None:
            raise ValueError("No sets were provided so the set size property was never set")
        total_pluckers = comb(num_cols,self.set_size)
        if self.which_positive is None:
            num_nonzeros = total_pluckers - len(self.which_zero)
        else:
            num_nonzeros = len(self.which_positive)
        if num_nonzeros == 0:
            raise ValueError(" ".join(["The constraints should be like a matroid specifying",
                                       "some (but not all) plucker coordinates as zero",
                                       "and the rest as nonzero"]))
        if sorter is None:
            # pylint: disable=unnecessary-lambda-assignment
            sorter = lambda z : z
        def my_set_leq(set_a : Set[int], set_b : Set[int]) -> bool:
            """
            set_a sorted compared to set_b sorted as lists
            """
            list_a = list(set_a)
            list_a.sort(key = sorter)
            list_b = list(set_b)
            list_b.sort(key = sorter)
            return list_a <= list_b
        if self.which_positive is not None:
            lex_min_set = reduce(lambda acc,x: x.copy() if my_set_leq(x,acc) else acc,
                                self.which_positive)
            return lex_min_set
        if self.which_zero is not None:
            which_positive = itertools.filterfalse(lambda item: item in self.which_zero,
                                map(set,
                                    itertools.combinations(range(self.n_at_least),self.set_size)
                                ))
            lex_min_set = reduce(lambda acc,x: x.copy() if my_set_leq(x,acc) else acc,
                                which_positive)
            return lex_min_set
        raise ValueError("At least one of which_positive and which_zero must be known")

    def to_hypersimplex_positroid_polytope(self,num_cols : int) -> Optional[ConvexHull]:
        """
        give the positroid polytope in the hypersimplex
        given by the image of the moment map from this positroid in tnn grassmanian
        to the hypersimplex
        See https://arxiv.org/abs/2104.08254
        """
        #pylint:disable=unused-variable
        def indicator_vector(which_is : Set[int]) -> np.ndarray:
            to_return = np.zeros(num_cols)
            for cur_i in which_is:
                to_return[cur_i] = 1
            return to_return
        all_pts = (indicator_vector(z) for z in self.which_positive)
        return create_convex_hull(all_pts)

    def cyclic_shift(self,my_n : int,num_times : int) -> None:
        """
        i->i-1
        0->n-1
        in all the sets
        """
        if self.which_positive is not None:
            for idx_in_list,cur_set in enumerate(self.which_positive):
                self.which_positive[idx_in_list] = {(z-num_times)%my_n for z in cur_set}
        if self.which_zero is not None:
            for idx_in_list,cur_set in enumerate(self.which_zero):
                self.which_zero[idx_in_list] = {(z-num_times)%my_n for z in cur_set}

    def check(self,k_subset : Set[int], cur_expr : Expr, var_set : Set[Symbol]) -> bool:
        """
        whichever constraint k_subset is in, check that sign constraint is satisfied
        """
        if self.which_positive is not None:
            if k_subset in self.which_positive:
                return nn_with_my_symbols(cur_expr,True,var_set)
            return cur_expr == Integer(0)
        if self.which_zero is not None:
            if k_subset in self.which_zero:
                return cur_expr == Integer(0)
            return nn_with_my_symbols(cur_expr,True,var_set)
        raise ValueError("At least one of which_positive and which_zero must be known")

class TNNGrassChart():
    """
    An image of R^{d}_{gt 0} in Mat(k,n)
    the Plucker coordinates have some sign constraints
    the purpose of this is for positroid charts
    """
    #pylint:disable=too-many-arguments
    def __init__(self, my_matrix : List[List[Expr]],*,variables : List[Symbol],
                 sign_constraints : Union[MinorSignConstraints,PositroidySignConstraints],
                 check_tnn : bool = True,gives_positroid_cell : bool = False):
        """
        the k by n (k<=n) matrix is given symbolically with d symbol variables
        listed in variables
        we can specify constraints on the maximal minors of the matrix
            the constraints are specified by which columns to pick out
        check_tnn specifies if we should explicitly make sure
            the sign constraints are satisfied in my_matrix
        """
        self.my_k = len(my_matrix)
        if self.my_k == 0:
            raise ValueError(f"The number of rows {self.my_k} must be positive")
        self.my_n = len(my_matrix[0])
        if self.my_k >= self.my_n:
            #pylint:disable=line-too-long
            raise ValueError(f"The number of rows {self.my_k} must be less than the number of columns {self.my_n}")
        if len(variables) > self.my_k*(self.my_n-self.my_k):
            raise ValueError("Dimension of this stratum is too big")
        self.variables = variables
        self.my_matrix = deepcopy(my_matrix)
        if sign_constraints.n_at_least>=self.my_n:
            raise ValueError("There were sign constraints with columns that do not exist")
        if sign_constraints.set_size != self.my_k:
            raise ValueError(f"Each constraint should have picked out {self.my_k} columns")
        if sign_constraints.is_positroidy(self.my_n):
            self.gives_positroid_cell = gives_positroid_cell
            if isinstance(sign_constraints,MinorSignConstraints):
                sign_constraints = sign_constraints.to_manifestly_positroidy(self.my_n)
        else:
            raise ValueError("The sign constraints were not the sort used for positroid cells")
        self.sign_constraints = sign_constraints
        self.cached_pluckers : Dict[FrozenSet,Expr] = {}
        if check_tnn:
            self.validate_pluckers()

    def cyclic_shift(self, num_times : int = 1) -> None:
        """
        on the matrix move the columns over by 1 with the first column becoming
        the last column after a sign correction with (-1)^(k-1)
        """
        if self.my_k % 2 == 0:
            num_times = num_times % self.my_n
        else:
            num_times = num_times % (2*self.my_n)
        if num_times == 0:
            return
        new_matrix = [[self.my_matrix[idx][(jdx+num_times) % self.my_n]
                       for jdx in range(self.my_n)]
                       for idx in range(self.my_k)]
        self.my_matrix = new_matrix
        if self.my_k % 2 == 0:
            for cur_col_from_end in range(1,num_times+1):
                for idx in range(self.my_k):
                    self.my_matrix[idx][(-cur_col_from_end) % self.my_n] *= Integer(-1)
        self.sign_constraints.cyclic_shift(self.my_n,num_times)
        new_dict : Dict[FrozenSet,Expr] = {}
        for key,value in self.cached_pluckers.items():
            new_key = frozenset((z-num_times)%self.my_n for z in key)
            new_dict[new_key] = value
        self.cached_pluckers = new_dict

    def validate_pluckers(self) -> None:
        """
        evaluate all the Plucker coordinates
        the ones in which_positive must be positive
        the rest must simply be nonnegative
        """
        var_set = set(self.variables)
        for k_subtuple in itertools.combinations(range(self.my_n),self.my_k):
            k_subset = set(k_subtuple)
            cur_plucker = self.get_plucker(k_subset)
            assert self.sign_constraints.check(k_subset,cur_plucker,var_set)

    def get_plucker(self,which_cols : Set[int]) -> Expr:
        """
        the Plucker coordinate for the specified k element subset of n columns
        """
        if len(which_cols) != self.my_k:
            raise ValueError(f"There should be {self.my_k} columns picked out")
        cached_val = self.cached_pluckers.get(frozenset(which_cols),None)
        if cached_val is not None:
            return cached_val
        which_cols_list = list(which_cols)
        which_cols_list.sort()
        my_matrix : List[List[Expr]] = [
            [self.my_matrix[row_idx][col_idx] for row_idx in range(self.my_k)]
                for col_idx in which_cols_list]
        answer = determinant(my_matrix)
        self.cached_pluckers[frozenset(which_cols)] = answer
        return answer

    def includes_in_which_schubert(self) -> Optional[Set[int]]:
        """
        which Schubert cell does this chart include into
        """
        try:
            to_return = self.sign_constraints.which_schubert_cell(self.my_n)
        except ValueError:
            to_return = None
        return to_return
