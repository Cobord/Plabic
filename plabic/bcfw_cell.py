"""
BCFW cells
https://arxiv.org/pdf/2402.15568.pdf
"""

from __future__ import annotations
from typing import Set

from plabic.plabic_diagram import PlabicGraph
from plabic.tnn_grassmannian import PositroidySignConstraints

#pylint:disable=too-few-public-methods
class ButterflyData:
    """
    the data needed to do a buttefly product
    """
    #pylint:disable=too-many-arguments
    def __init__(self,a_gluing : int,b_gluing : int,
                 c_gluing : int,d_gluing : int,n_gluing : int):
        self.a_gluing = a_gluing
        self.b_gluing = b_gluing
        self.c_gluing = c_gluing
        self.d_gluing = d_gluing
        self.n_gluing = n_gluing

class BCFWCell:
    """
    BCFW cells are recursively defined
    - the cells represented by
        the PlabicGraph instances which only have
        Color.?? (TODO check which one I called Red vs Green) lollipops
        at all of the boundary vertices
    - rotations, reflections and zero insertions of some other BCFW cell
    - a butterfly product along a ButteflyData 
    """
    #pylint:disable=unused-private-member
    __my_kaleidoscope : PlabicGraph # because collection of butterflies = kaleidoscope
    __my_positroid : PositroidySignConstraints
    __is_manifestly_standard : bool

    def __init__(self, k : int):
        """
        the base cases
        """
        raise NotImplementedError

    def rotate(self, how_many_times : int = 1) -> None:
        """
        mutate self to a rotation of itself
        """
        if how_many_times == 0:
            return
        self.__is_manifestly_standard = False
        raise NotImplementedError

    def reflect(self) -> None:
        """
        mutate self to a reflection of itself
        """
        self.__is_manifestly_standard = False
        raise NotImplementedError

    def zero_insertions(self, j_set : Set[int]) -> None:
        """
        insert zero columns
        """
        if j_set.size() == 0:
            return
        raise NotImplementedError

    def butterfly_product(self,other : BCFWCell, _along_this : ButterflyData) -> None:
        """
        replace self with self X_along_this other
        """
        #pylint:disable=protected-access
        self.__is_manifestly_standard = self.__is_manifestly_standard and \
            other.__is_manifestly_standard
        raise NotImplementedError
