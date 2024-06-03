"""
BCFW cells
https://arxiv.org/pdf/2402.15568.pdf
"""

#pylint:disable=line-too-long,too-many-locals,unused-variable,unreachable

from __future__ import annotations
import random
import string
from typing import List, Set

from .plabic_diagram import PlabicGraph, PlabicGraphBuilder
from .tnn_grassmannian import PositroidySignConstraints

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
    __my_kaleidoscopes_outer_based_labelling : List[str]
    __my_positroid : PositroidySignConstraints
    __is_manifestly_standard : bool

    def __init__(self, k : int):
        """
        the base cases
        first bullet point above
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

    def butterfly_product(self,other : BCFWCell, along_this : ButterflyData) -> None:
        """
        replace self with self X_along_this other
        """
        #pylint:disable=protected-access
        self.__is_manifestly_standard = self.__is_manifestly_standard and \
            other.__is_manifestly_standard
        butterfly_template_builder = PlabicGraphBuilder()
        gl_externals = self.__my_kaleidoscopes_outer_based_labelling
        gr_externals = other.__my_kaleidoscopes_outer_based_labelling
        butterfly_template_builder.set_internal_circles_nums([len(gl_externals),len(gr_externals)])
        butterfly_template_builder.set_num_external(len(gl_externals)+len(gr_externals)-1)
        my_based_labelling = ["" for _ in range(len(gl_externals)+len(gr_externals)-1)]
        cur_external_bdry_idx = 0
        # the lines connecting 1 through a-1 directly into G_L
        shift_amt = len(gl_externals) - along_this.b_gluing - 1
        for bdry_idx in range(1,along_this.a_gluing):
            name_in_gl = gl_externals[(bdry_idx - 1 - shift_amt) % len(gl_externals) ]
            new_name = name_in_gl + "".join(random.sample(string.ascii_letters,5))
            butterfly_template_builder.add_external_bdry_vertex(new_name,cur_external_bdry_idx,name_in_gl)
            my_based_labelling[cur_external_bdry_idx] = new_name
            cur_external_bdry_idx += 1
        gl_a_name = gl_externals[(along_this.a_gluing - 1 - shift_amt) % len(gl_externals) ]
        gl_b_name = gl_externals[(along_this.a_gluing - shift_amt) % len(gl_externals) ]
        gl_n_name = gl_externals[(along_this.a_gluing + 1 - shift_amt) % len(gl_externals) ]
        # the lines connecting b+1 through c-1 directly into G_R
        cur_external_bdry_idx += 2
        bs_idx_in_gr = (along_this.n_gluing + 1) % len(gr_externals)
        for bdry_idx in range(along_this.b_gluing+1,along_this.c_gluing):
            name_in_gr = gr_externals[(bdry_idx - along_this.b_gluing + \
                                       bs_idx_in_gr) % len(gr_externals) ]
            new_name = name_in_gr + "".join(random.sample(string.ascii_letters,5))
            butterfly_template_builder.add_external_bdry_vertex(new_name,cur_external_bdry_idx,name_in_gl)
            my_based_labelling[cur_external_bdry_idx] = new_name
            cur_external_bdry_idx += 1
        gr_c_name = gr_externals[(along_this.c_gluing - along_this.b_gluing + \
                                  bs_idx_in_gr) % len(gr_externals) ]
        gr_d_name = gr_externals[(along_this.c_gluing - along_this.b_gluing + 1 + \
                                  bs_idx_in_gr) % len(gr_externals) ]
        gr_n_name = gr_externals[(along_this.c_gluing - along_this.b_gluing + 2 + \
                                  bs_idx_in_gr) % len(gr_externals) ]
        gr_b_name = gr_externals[(along_this.c_gluing - along_this.b_gluing + 3 + \
                                  bs_idx_in_gr) % len(gr_externals) ]
        raise NotImplementedError
        butterfly = butterfly_template_builder.build()
        success, reason = butterfly.operad_compose(self.__my_kaleidoscope,0)
        if not success:
            return reason
        success, reason = butterfly.operad_compose(other.__my_kaleidoscope,0)
        if not success:
            return reason
        self.__my_kaleidoscope = butterfly
        self.__my_kaleidoscopes_outer_based_labelling = my_based_labelling
        # self.__my_positroid = ???
