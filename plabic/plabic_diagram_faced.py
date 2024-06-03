"""
PlabicGraph class does not store the faces
this is used on top of that
mostly delegates mutation stuff to the underlying PlabicGraph
completely delegates the query stuff to the underlying PlabicGraph
"""

from __future__ import annotations
from itertools import chain
from typing import Callable, Iterable, List, Optional, Tuple, TypeVar

from sympy import Expr, Integer, symbols
import networkx as nx

from .cluster import Arithmetic, Cluster
from .cyclic_utils import combine_cyclicly_ordered, put_in_between
from .plabic_diagram import BiColor, ExtraData, PlabicGraph, Point

NT = TypeVar("NT",bound=Arithmetic)

#pylint:disable=too-many-public-methods
class PlabicDiagramFaced:
    """
    PlabicGraph avoids storing the faces because that does not matter for rules of the road/drawing
    so separate concerns to this class that encapsulates PlabicGraph
    """
    def __init__(self, pb: PlabicGraph, faces: Iterable[Iterable[str]]):
        self.__pb = pb
        self.__faces = [list(face) for face in faces]
        try:
            initial_variables : List[Expr] = list(symbols(",".join(
                [f"x_{i}" for i in range(len(self.__faces))]),
                commutative=True))
            _junk : Cluster[Expr] = self.to_cluster(initial_variables,
                   Integer(1),
                   lambda z: z.simplify())
        except ValueError:
            #pylint:disable=raise-missing-from
            raise ValueError(
                "Something wrong with faces data provided as describing the faces in pb")

    def is_trivalent(self, this_vertex: str) -> bool:
        """
        is this_vertex trivalent
        """
        return self.__pb.is_trivalent(this_vertex)

    def one_color_trivalent(self, all_of_this_color: BiColor) -> bool:
        """
        are all vertices of specified color trivalent
        """
        return self.__pb.one_color_trivalent(all_of_this_color)

    def boundary_connects_to_color(self, all_of_this_color: BiColor) -> bool:
        """
        are boundary vertices connected to only internal vertices of specified color
        no boundary to boundary allowed, no boundary to the opposite color allowed
        """
        return self.__pb.boundary_connects_to_color(all_of_this_color)

    def is_bipartite(self) -> bool:
        """
        are all edges connecting internal vertices opposite colors
        """
        return self.__pb.is_bipartite()

    def nodes_connected(self, this_vertex: str, that_vertex: str) -> bool:
        """
        are this_vertex and that_vertex connected by an edge
        """
        return self.__pb.nodes_connected(this_vertex, that_vertex)

    def num_connecting_edges(self, this_vertex: str, that_vertex: str) -> int:
        """
        how many edges connect this_vertex and that_vertex
        """
        return self.__pb.num_connecting_edges(this_vertex, that_vertex)

    def get_color(self, this_vertex: str) -> Optional[BiColor]:
        """
        the color of this vertex
        """
        return self.__pb.get_color(this_vertex)

    def opposite_colors(self, this_vertex: str, that_vertex: str) -> bool:
        """
        are this_vertex and that_vertex both colored and opposite
        """
        return self.__pb.opposite_colors(this_vertex, that_vertex)

    def same_colors(self, this_vertex: str, that_vertex: str) -> bool:
        """
        are this_vertex and that_vertex both colored and the same
        """
        return self.__pb.same_colors(this_vertex, that_vertex)

    def my_kn_type(self, this_contribution : Optional[str]=None) -> Tuple[int,int]:
        """
        bdry_sinks - bdry_sources + 
            sum_{v in red vertices} (deg(v)-2)
            - sum_{v in green vertices} (deg(v)-2) = 0
        bdry_sinks - bdry_sources = - my_sum
        bdry_sinks + bdry_sources = my_n
        """
        return self.__pb.my_kn_type(this_contribution)

    def greedy_shrink(self,
                      name_combiner : Optional[Callable[[str,str],str]] = None,
                      rounds_left : int = 5) -> bool:
        """
        similar to greedy_shrink of underlying PlabicGraph
        but using own versions of contract_edge and remove_bivalent
        so that faces are taken care of too
        """
        if not self.__pb.my_extra_props.issubset([]):
            return False
        name_combiner = name_combiner if name_combiner is not None else lambda z1,_: z1
        any_change = False
        all_edges = list(self.__pb.my_graph.edges(keys=False))
        for name1,name2 in all_edges:
            if name1 in self.__pb.my_graph.nodes() and name2 in self.__pb.my_graph.nodes():
                did_remove,_ = self.contract_edge(name1,name2,name_combiner(name1,name2))
            else:
                did_remove = False
            any_change = any_change or did_remove
        all_node_names = list(self.__pb.my_graph.nodes())
        for name in all_node_names:
            did_remove,_ = self.remove_bivalent(name)
            any_change = any_change or did_remove
        if not any_change:
            return False
        if rounds_left>0:
            _ = self.greedy_shrink(name_combiner,rounds_left-1)
        return True

    def square_move(self, four_nodes: Tuple[str, str, str, str]) -> Tuple[bool, str]:
        """
        purely delegate the square move
        no need to change __faces because names of vertices don't change
        """
        return self.__pb.square_move(four_nodes)

    # pylint:disable = too-many-return-statements, too-many-statements
    def flip_move(self, this_vertex: str,
                  that_vertex: str,
                  extra_data_transformer:
                  Optional[Callable[[ExtraData, ExtraData, "PlabicGraph"],
                                    Tuple[ExtraData, ExtraData]]] = None)\
            -> Tuple[bool, str]:
        """
        requires implementation details of PlabicGraph's flip_move
        """
        raise NotImplementedError

    def remove_bivalent(self, my_bivalent_vertex: str) -> Tuple[bool, str]:
        """
        remove the specified bivalent vertex
        if these names did not meet the pre-conditions of being bivalent
            provide the explanation as well
        also output if self has changed
        """
        success, explanation = self.__pb.remove_bivalent(my_bivalent_vertex)
        if not success:
            return success, explanation
        for face in self.__faces:
            try:
                face.remove(my_bivalent_vertex)
            except ValueError:
                pass
        return success, explanation

    #pylint:disable=too-many-arguments
    def insert_bivalent(self, this_vertex: str, that_vertex: str,
                        desired_color: BiColor,
                        desired_name: str,
                        extra_data_transformer:
                            Optional[Callable[[ExtraData, ExtraData, "PlabicGraph"],
                                              ExtraData]] = None) -> Tuple[bool, str]:
        """
        add a bivalent vertex on the edge connecting these two vertices
        if these names did not meet the pre-conditions of having an edge
            provide the explanation as well
        also output if self has changed
        """
        for idx,face in enumerate(self.__faces):
            if this_vertex in face and that_vertex in face:
                this_idx = face.index(this_vertex)
                that_idx = face.index(that_vertex)
                if abs(this_idx-that_idx) == 1:
                    continue
                if this_idx==0 and that_idx==face.len()-1:
                    continue
                if that_idx==0 and this_idx==face.len()-1:
                    continue
                return False, " ".join(["The two vertices should be connected",
                                       "by a single edge so any face that has both,",
                                       "should have them next to each other"])
        success, explanation = \
            self.__pb.insert_bivalent(this_vertex,that_vertex,
                                      desired_color,desired_name,extra_data_transformer)
        if not success:
            return success, explanation
        for idx,face in enumerate(self.__faces):
            if this_vertex in face and that_vertex in face:
                self.__faces[idx] = put_in_between(face, this_vertex, that_vertex, desired_name)
        return success, explanation

    def contract_edge(self, this_vertex: str,
                      that_vertex: str,
                      combined_name: str,
                      extra_data_transformer:
                      Optional[Callable[[ExtraData, ExtraData, "PlabicGraph"],
                                        ExtraData]] = None)\
            -> Tuple[bool, str]:
        """
        contract the edge connecting these two vertices of the same color
        if these names did not meet the pre-conditions
            provide the explanation as well
        also output if self has changed
        """
        for idx,face in enumerate(self.__faces):
            if this_vertex in face and that_vertex in face:
                this_idx = face.index(this_vertex)
                that_idx = face.index(that_vertex)
                if abs(this_idx-that_idx) == 1:
                    continue
                if this_idx==0 and that_idx==face.len()-1:
                    continue
                if that_idx==0 and this_idx==face.len()-1:
                    continue
                return False, " ".join(["The two vertices should be connected by",
                                        "a single edge so any face that has both,",
                                        "should have them next to each other"])
        success, explanation = \
            self.__pb.contract_edge(this_vertex,that_vertex,
                                    combined_name,extra_data_transformer)
        if not success:
            return success, explanation
        for idx, face in enumerate(self.__faces):
            if this_vertex in face and that_vertex in face:
                face.remove(this_vertex)
                where_that = face.index(that_vertex)
                face[where_that] = combined_name
                self.__faces[idx] = face
            elif this_vertex in face:
                where_this = face.index(this_vertex)
                face[where_this] = combined_name
                self.__faces[idx] = face
            elif that_vertex in face:
                where_that = face.index(that_vertex)
                face[where_that] = combined_name
                self.__faces[idx] = face
        return success, explanation

    #pylint:disable=too-many-arguments
    def split_vertex(self, this_vertex: str,
                     split_bounds: Tuple[int, int],
                     split_name_1: str,
                     split_name_2: str,
                     extra_data_transformer:
                     Optional[Callable[[ExtraData, "PlabicGraph"],
                                       Tuple[ExtraData, ExtraData]]] = None)\
            -> Tuple[bool, str]:
        """
        split this vertex into an edge connecting two vertices of the same color
        split_bounds informs which edges attach to which
            it splits the cyclic ordering into two linear orderings
        if anything did not meet the pre-conditions
            provide the explanation as well
        also output if self has changed
        """
        raise NotImplementedError

    def radius_unoccupied(self, center : Point, ignored_vertices : List[str]) -> Tuple[bool,float]:
        """
        the ball of the returned radius around center
        is only allowed to have ignored_vertices within it
        """
        return self.__pb.radius_unoccupied(center, ignored_vertices)

    def coordinate_transform(self,transform : Callable[[Point],Point],
                             circle_config_invalidates : bool = True) -> bool:
        """
        does transform on all positions of all nodes
        """
        return self.__pb.coordinate_transform(transform,circle_config_invalidates)

    #pylint:disable=too-many-locals
    def operad_compose(self,
                       other: PlabicDiagramFaced,
                       which_internal_disk: int) -> Tuple[bool, str]:
        """
        substitute other in on the i'th internal circle of self
            that is self has at least that many in self.my_internal_bdry
            and self.my_internal_bdry[which_internal_disk] are matched up
            with the other.my_external_nodes
        if anything did not meet the pre-conditions
            provide the explanation as well
        also output if self has changed
        """
        glued_vertices = self.__pb.my_internal_bdry[which_internal_disk]
        unaffected_self_faces = [face for face in self.__faces
                                 if all(g not in face for g in glued_vertices)]
        #pylint:disable=protected-access
        unaffected_other_faces = [face for face in other.__faces
                                  if all(g not in face for g in glued_vertices)]
        to_iterate = zip(glued_vertices,glued_vertices[1:])
        to_iterate = chain(to_iterate, [(glued_vertices[-1], glued_vertices[0])])
        cross_faces : List[List[str]] = []
        for glued_vertex_1,glued_vertex_2 in to_iterate:
            from_self = [face for face in self.__faces
                         if glued_vertex_1 in face and glued_vertex_2 in face]
            from_other = [face for face in other.__faces
                          if glued_vertex_1 in face and glued_vertex_2 in face]
            success, new_face, message = \
                combine_cyclicly_ordered(from_self,from_other, glued_vertex_1, glued_vertex_2)
            if not success:
                return False, message
            cross_faces.append(new_face)
        change_made, reason = self.__pb.operad_compose(other.__pb, which_internal_disk)
        if not change_made:
            return change_made, reason
        self.__faces = list(chain(unaffected_other_faces,unaffected_self_faces,cross_faces))
        return change_made, reason

    def bdry_to_bdry(self, starting_bdry_node: str,
                     care_about_lollipop : bool = True) -> Tuple[Optional[BiColor], str]:
        """
        starting at a boundary vertex follow the rules of the road until another
        boundary vertex (possibly the same one)
        the color gives the decorated permutation
        because fixed points are colored by the color of the lollipop
        """
        return self.__pb.bdry_to_bdry(starting_bdry_node, care_about_lollipop)

    def follow_rules_of_road(self, at_this_node: str,
                             was_at_this_before: Optional[str] = None,
                             via_edge_keyed: Optional[int] = None) -> Optional[Tuple[int, str]]:
        """
        if we are at some vertex and came from some other node, where should we go next
        if we are following the rule that we turn maximally right/left
        depending on the color of where we are and stop if we are on a boundary vertex
        """
        return self.follow_rules_of_road(at_this_node, was_at_this_before, via_edge_keyed)

    def draw(self, *,
             draw_oriented_if_perfect=True,
             show_node_names : bool = True,
             red_nodes: str = "red",
             green_nodes: str = "green",
             bdry_nodes: str = "black",
             oriented_arrows: str = "black",
             unoriented_arrows_perfect: str = "yellow",
             unoriented_arrows_imperfect: str = "black",
             overridden_arrow_orientation : Optional[Callable[[str,str,int],bool]] = None,
             external_circle_color : str = "black",
             internal_circles_color : str = "black",
             show_as_well : bool = True
             ) -> None:
        """
        draw the multigraph without regard to it's planar embedding

        draw_oriented_if_perfect - if the Plabic graph is equipped with a perfect
            orientation, draw the edges as arrows instead of line segments
            if the graph is not bipartite, then drawing as line segments does not work
            so perfect orientations will be drawn with arrows regardless of this input
        show_node_names - show the names of the vertices
        red_nodes,green_nodes,bdry_nodes - the colors of the vertices
            (adjust for black and white printing or red green colorblindness)
        oriented_arrows - what color to draw oriented edges
        unoriented_arrows_perfect,unoriented_arrows_imperfect - if the Plabic graph
            is equipped with a perfect orientation and we are drawing edges without
            orientation then what color to mark the perfect/imperfect edges
        overridden_arrow_orientation - a function which takes the two endpoints
            and key (for multiedges) and if the output is false
            then it doesn't draw that arrowhead for that orientation of the edge
        external_circle_color,internal_circles_color - if the boundaries are on circles
            in a FramedDiskConfig then what color to draw those framed circles
        show_as_well - show right away or let the caller do the plt.show
        """
        return self.__pb.draw(draw_oriented_if_perfect=draw_oriented_if_perfect,
             show_node_names=show_node_names,
             red_nodes=red_nodes,
             green_nodes=green_nodes,
             bdry_nodes=bdry_nodes,
             oriented_arrows=oriented_arrows,
             unoriented_arrows_perfect=unoriented_arrows_perfect,
             unoriented_arrows_imperfect=unoriented_arrows_imperfect,
             overridden_arrow_orientation=overridden_arrow_orientation,
             external_circle_color=external_circle_color,
             internal_circles_color=internal_circles_color,
             show_as_well=show_as_well)

    def to_cluster(self,
                   initial_variables: List[NT],
                   multiplicative_identity: NT,
                   simplifier: Callable[[NT],NT]) -> Cluster[NT]:
        """
        construct quiver with vertices specifed by self's faces
        """
        count_vertices = len(self.__faces)
        my_quiver = nx.MultiDiGraph()
        for cur_node in range(1,count_vertices+1):
            my_quiver.add_node(str(cur_node))
        for this_face in range(1,count_vertices+1):
            this_face_vertices = self.__faces[this_face]
            for that_face in range(this_face+1,count_vertices+1):
                that_face_vertices = self.__faces[that_face]
                if any(t in this_face_vertices for t in that_face_vertices):
                    # add edges
                    raise NotImplementedError
        return Cluster[NT](my_quiver,initial_variables, multiplicative_identity, simplifier)
