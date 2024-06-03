"""
plabic graph
PLAnar BIColored
"""
# pylint:disable=too-many-lines
from __future__ import annotations
from enum import Enum, auto
import functools
import numbers
from math import sqrt
from typing import Tuple, Optional, List, Dict, cast, Any, Callable, Set, Iterator
import itertools
import networkx as nx
import matplotlib.pyplot as plt

from .cyclic_utils import cyclic_equal
from .framed_2_disks import FramedDiskConfig

Point = Tuple[float,float]

class BiColor(Enum):
    """
    the two colors for internal vertices
    """
    RED = auto()
    GREEN = auto()


ExtraData = Dict[str, Any]

# pylint:disable = too-many-instance-attributes
class PlabicGraphBuilder:
    """
    PlabicGraph constructor arguments are complicated to offload this process
    to PlabicGraphBuilder
    """
    my_init_data : Dict[str, Tuple[BiColor, List[str]]]
    num_external_vertices : int
    is_set_externals : List[bool]
    external_init_orientation : List[str]
    multi_edge_permutation : Dict[Tuple[str, str], Dict[int, int]]
    internal_circles_nums : List[int]
    is_set_internals : List[List[bool]]
    internal_bdry_orientations: Optional[List[List[str]]]
    extra_node_props: Optional[Dict[str, ExtraData]]
    have_multi_edges : Set[Tuple[str,str]]
    circles_config : Optional[FramedDiskConfig]
    all_names_anywhere : Set[str]

    def __init__(self):
        self.my_init_data = {}
        self.num_external_vertices = 0
        self.is_set_externals = []
        self.external_init_orientation = []
        self.multi_edge_permutation = {}
        self.internal_circles_nums = []
        self.is_set_internals = []
        self.internal_bdry_orientations = None
        self.extra_node_props = None
        self.have_multi_edges = set()
        self.circles_config = None
        self.all_names_anywhere = set()

    def set_num_external(self,num_ext : int) -> None:
        """
        set the number of vertices on the outside boundary
        """
        if num_ext<0:
            raise ValueError("Num external vertices should be nonnegative")
        self.num_external_vertices = num_ext
        self.external_init_orientation = [""]*num_ext
        self.is_set_externals = [False]*num_ext

    def set_internal_circles_nums(self,internal_circles_nums : List[int]):
        """
        set the number of internal boundaries and how many vertices are on each
        """
        if any(z<0 for z in internal_circles_nums):
            raise ValueError("Num vertices on an internal boundary should all be nonnegative")
        self.internal_circles_nums = internal_circles_nums
        self.is_set_internals = [[False]*z for z in internal_circles_nums]
        self.internal_bdry_orientations = [[""]*z for z in internal_circles_nums]

    def set_extras(self, vertex_name : str, extras : Optional[ExtraData] = None):
        """
        set the extra_node_props for a specific vertex
        """
        self.all_names_anywhere.add(vertex_name)
        if extras is not None:
            if self.extra_node_props is None:
                self.extra_node_props = {}
            self.extra_node_props[vertex_name] = extras

    def add_internal_vertex(self,vertex_name : str, color : BiColor,
                            clockwise_nhbrs : List[str],
                            extras : Optional[ExtraData] = None) -> None:
        """
        add an internal vertex by it's color and it's half-edges around it
        the half-edges are in clockwise order, they give the vertex at the other
        end
        any extra information want on that node is also provided
        """
        self.all_names_anywhere.add(vertex_name)
        for name in clockwise_nhbrs:
            self.all_names_anywhere.add(name)
            if clockwise_nhbrs.count(name)>1:
                self.have_multi_edges.add((vertex_name,name))
                self.have_multi_edges.add((name,vertex_name))
        self.my_init_data[vertex_name] = (color,clockwise_nhbrs)
        self.set_extras(vertex_name,extras)

    def add_external_bdry_vertex(self,external_vertex_name : str,
                            idx : int, my_internal_connection : str,
                            extras : Optional[ExtraData] = None) -> None:
        """
        add a boundary vertex on the external boundary by it's only neighbor
        any extra information want on that node is also provided
        it is placed on the external boundary in position provided by idx
        you must have done set_num_external before this so we know what
        range idx is regarded as taking values in
        """
        if idx<0 or idx>=self.num_external_vertices:
            raise ValueError(f"{self.num_external_vertices} external vertices allowed")
        self.all_names_anywhere.add(external_vertex_name)
        self.all_names_anywhere.add(my_internal_connection)
        self.my_init_data[external_vertex_name] = (BiColor.RED,[my_internal_connection])
        self.external_init_orientation[idx] = external_vertex_name
        self.is_set_externals[idx] = True
        self.set_extras(external_vertex_name,extras)

    # pylint:disable=too-many-arguments
    def add_internal_bdry_vertex(self,internal_bdry_vertex_name : str,
                            which_circle : int,
                            where_in_that_circle : int,
                            my_internal_connection : str,
                            extras : Optional[ExtraData] = None) -> None:
        """
        add a boundary vertex on an internal boundary by it's only neighbor
        any extra information want on that node is also provided
        it is placed on the internal boundary which_circle in position provided by
        where_in_that_circle
        you must have done set_internal_circles_nums before this so we know what
        range which_circle and where_in_that_circle are regarded as taking values in
        """
        if self.internal_bdry_orientations is None:
            raise ValueError(
                "Should have set the number and sizes of the internal boundaries already")
        if which_circle<0 or which_circle>=len(self.internal_circles_nums):
            raise ValueError(
                f"{which_circle} is not one of the available internal circles")
        if where_in_that_circle<0 or where_in_that_circle>=self.internal_circles_nums[which_circle]:
            raise ValueError(
                f"{where_in_that_circle} is not available on internal circle {which_circle}")
        self.all_names_anywhere.add(internal_bdry_vertex_name)
        self.all_names_anywhere.add(my_internal_connection)
        self.my_init_data[internal_bdry_vertex_name] = (BiColor.RED,[my_internal_connection])
        self.internal_bdry_orientations[which_circle][where_in_that_circle] = \
            internal_bdry_vertex_name
        self.is_set_internals[which_circle][where_in_that_circle] = True
        self.set_extras(internal_bdry_vertex_name,extras)

    def set_multi_edges(self,vertex_a : str, vertex_b: str, how_keys_match : Dict[int,int]) -> None:
        """
        set the permutation of half-edges if there are multiple edges connecting
        vertex_a and vertex_b
        """
        self.all_names_anywhere.add(vertex_a)
        self.all_names_anywhere.add(vertex_b)
        self.multi_edge_permutation[(vertex_a,vertex_b)] = how_keys_match

    def set_circles_config(self,circles_config : FramedDiskConfig) -> None:
        """
        if all the positions of the boundary vertices are set,
        then we can set the circles that they all lie on
        """
        if len(circles_config.internal_circles) != len(self.internal_circles_nums):
            raise ValueError("There are not the right number of internal circles")
        if not all(self.is_set_externals):
            raise ValueError("Not finished setting externals")
        if self.is_set_internals is not None and not all(all(z) for z in self.is_set_internals):
            raise ValueError("Not finished setting internals")
        self.circles_config = circles_config

    # pylint:disable = too-many-locals
    def rename(self, old_name : str , new_name : str) -> None:
        """
        Change all occurences of vertex name old_name to new_name
        """
        if old_name not in self.all_names_anywhere:
            return
        if new_name in self.all_names_anywhere:
            raise ValueError(" ".join([
                f"{new_name} was already a vertex name, so this",
                "renaming would cause name collisions"]))

        def replace_in_list(name_list : List[str]) -> bool:
            """
            in place replacement in a list of old_name with new_name
            """
            any_changes = False
            for idx, item in enumerate(name_list):
                if item == old_name:
                    name_list[idx] = new_name
                    any_changes = True
            return any_changes
        def fix_name(a_name : str) -> str:
            """
            the modified name
            """
            return new_name if a_name==old_name else a_name

        try:
            old_data = self.my_init_data.pop(old_name)
            _ = replace_in_list(old_data[1])
            self.my_init_data[new_name] = (old_data[0],old_data[1])
        except KeyError:
            pass
        for _,value in self.my_init_data.items():
            _ = replace_in_list(value[1])
        _ = replace_in_list(self.external_init_orientation)
        needs_renamings = [k for k in self.multi_edge_permutation if old_name in k]
        for key in needs_renamings:
            old_value = self.multi_edge_permutation.pop(key)
            new_key = (fix_name(key[0]),fix_name(key[1]))
            self.multi_edge_permutation[new_key] = old_value
        if self.internal_bdry_orientations is not None:
            for cur_list in self.internal_bdry_orientations:
                _ = replace_in_list(cur_list)
        if self.extra_node_props is not None:
            try:
                old_extra_data = self.extra_node_props.pop(old_name)
                self.extra_node_props[new_name] = old_extra_data
            except KeyError:
                pass
        bad_multi_tuples = [k for k in self.have_multi_edges if old_name in k]
        for bad_tuple in bad_multi_tuples:
            good_tuple = (fix_name(bad_tuple[0]),fix_name(bad_tuple[1]))
            self.have_multi_edges.remove(bad_tuple)
            self.have_multi_edges.add(good_tuple)
        self.all_names_anywhere.remove(old_name)
        self.all_names_anywhere.add(new_name)

    def build(self) -> PlabicGraph:
        """
        use all the instance attributes in arguments to
        PlabicGraph construction
        """
        if not all(self.is_set_externals):
            raise ValueError("Not finished setting externals")
        if self.is_set_internals is not None and not all(all(z) for z in self.is_set_internals):
            raise ValueError("Not finished setting internals")
        for (vertex_a,vertex_b) in self.have_multi_edges:
            ab_in_multi = (vertex_a,vertex_b) in self.multi_edge_permutation
            ba_in_multi = (vertex_b,vertex_a) in self.multi_edge_permutation
            if not ab_in_multi and not ba_in_multi:
                raise ValueError("Not finished setting for multi-edges")
        for vertex_name in self.all_names_anywhere:
            if vertex_name not in self.my_init_data:
                raise ValueError("Not finished setting neighbors of all vertices")
        return PlabicGraph(self.my_init_data,
                 self.external_init_orientation,
                 self.multi_edge_permutation,
                 self.internal_bdry_orientations,
                 self.extra_node_props,
                 self.circles_config)

# pylint:disable=too-many-public-methods,too-many-instance-attributes
class PlabicGraph:
    """
    built on top of a directed multigraph
    """

    num_interior_circles: int
    my_graph: nx.MultiDiGraph
    my_external_nodes: List[str]
    my_internal_bdry: List[List[str]]
    all_bdry_nodes: List[str]
    multi_edge_permutation: Dict[Tuple[str, str], Dict[int, int]]
    my_perfect_matching: Optional[Set[Tuple[str, str, int]]]
    my_extra_props : Set[str]
    circles_config : Optional[FramedDiskConfig]

    # pylint:disable = too-many-arguments, too-many-locals, too-many-branches,too-many-statements
    def __init__(self, my_init_data: Dict[str, Tuple[BiColor, List[str]]],
                 external_init_orientation: List[str],
                 multi_edge_permutation: Dict[Tuple[str, str], Dict[int, int]],
                 internal_bdry_orientations: Optional[List[List[str]]] = None,
                 extra_node_props: Optional[Dict[str, ExtraData]] = None,
                 my_circles_config : Optional[FramedDiskConfig] = None):
        """
        provide the data of the graph by a dictionary
        the keys are the names of the vertices
        the values are a tuple of it's color if it is internal
        and the clockwise list of it's neighboring nodes
        the list only matters up to rotation
        external_orientation is the list of boundary nodes and their indices
        give the 1 through n clockwise labels
        """
        if extra_node_props is None:
            extra_node_props = {}
        if internal_bdry_orientations is None:
            self.num_interior_circles = 0
            internal_bdry_orientations = []
        else:
            self.num_interior_circles = len(internal_bdry_orientations)
        self.my_graph = nx.MultiDiGraph()
        all_bdry_vertices = external_init_orientation.copy()
        for internal_circle_num in range(self.num_interior_circles):
            all_bdry_vertices.extend(
                internal_bdry_orientations[internal_circle_num])
        if len(set(all_bdry_vertices)) != len(all_bdry_vertices):
            raise ValueError("The boundary vertices must all be distinct")
        if any(z not in my_init_data for z in all_bdry_vertices):
            raise ValueError("The boundary vertices must all have data")
        self.my_extra_props : Set[str] = set()
        self.my_external_nodes = external_init_orientation
        self.my_internal_bdry = internal_bdry_orientations
        self.all_bdry_nodes = all_bdry_vertices
        for src, (cur_color, tgts) in my_init_data.items():
            is_cur_bdry = src in all_bdry_vertices
            if is_cur_bdry:
                self.my_graph.add_node(src, is_interior=False)
            else:
                self.my_graph.add_node(src, is_interior=True,
                                       color=cur_color)
            self.__add_props(src, extra_node_props.get(src, {}))
            if is_cur_bdry and len(tgts) != 1:
                raise ValueError(
                    "Each boundary vertex must be incident to a single edge")
            for edge_idx, tgt in enumerate(tgts):
                self.my_graph.add_edge(src, tgt, key=edge_idx)
        node_list = enumerate(my_init_data.keys())
        self.multi_edge_permutation: Dict[Tuple[str, str], Dict[int, int]] = {}
        any_self_loops = False
        for (src_num, src), (tgt_num, tgt) in itertools.combinations_with_replacement(node_list, 2):
            if src_num > tgt_num:
                continue
            if src_num == tgt_num:
                num_src_loops = self.num_connecting_edges(src, src)
                if num_src_loops == 0:
                    continue
                any_self_loops = True
                if num_src_loops > 1:
                    src_to_src_dict = multi_edge_permutation.get(
                        (src, src), None)
                    self.__multi_edge_permutation_add(
                        src, src, src_to_src_dict, None, num_src_loops)
                else:
                    msg = " ".join([f"If a self loop at {src},",
                                    "then there would be at least 2",
                                    "half edges there"])
                    raise ValueError(msg)
                continue
            num_src_tgt = self.num_connecting_edges(src, tgt)
            num_tgt_src = self.num_connecting_edges(tgt, src)
            if num_src_tgt != num_tgt_src:
                raise ValueError(
                    "\n".join([
                        "Each edge needs to have a corresponding backwards edge",
                        f"The one currently wrong was {src} {tgt}",
                        f"Forward {num_src_tgt}, Backwards {num_tgt_src}",
                        f"All out of {src} : {my_init_data[src][1]}",
                        f"All out of {tgt} : {my_init_data[tgt][1]}",
                        ]))
            if num_src_tgt > 1:
                src_to_tgt_dict = multi_edge_permutation.get((src, tgt), None)
                tgt_to_src_dict = multi_edge_permutation.get((tgt, src), None)
                self.__multi_edge_permutation_add(src, tgt, src_to_tgt_dict,
                                                  tgt_to_src_dict, num_src_tgt)
        self.__boundary_connectivity()
        self.my_perfect_matching: Optional[Set[Tuple[str, str, int]]] = set()
        self.__my_perfect_matching_fix(any_self_loops)

        self.circles_config = None
        if my_circles_config is not None:
            self.circles_config = my_circles_config
            # trusts that all the boundary vertices on the circles they are supposed to be

    def __add_props(self, vertex: str, props: ExtraData):
        """
        add all the extra keys present for extra data to attach to nodes
        for example, position
        """
        for prop_name, prop_val in props.items():
            if prop_name not in ["is_interior", "color"]:
                self.my_graph.nodes[vertex][prop_name] = prop_val
                self.my_extra_props.add(prop_name)

    def remove_prop(self, prop_name : str):
        """
        remove a property from all nodes
        """
        if prop_name in ["is_interior", "color"]:
            raise ValueError("Can only delete the optional properties")
        for cur_node in self.my_graph.nodes:
            try:
                del self.my_graph.nodes[cur_node][prop_name]
            except KeyError:
                pass
        self.my_extra_props.discard(prop_name)

    def __boundary_connectivity(self):
        """
        make sure every internal vertex is connected to some external boundary vertex
        """
        infinity_node = "infinity"
        assert infinity_node not in self.my_graph
        self.my_graph.add_node(
            infinity_node, is_interior=False, color=BiColor.RED)
        for edge_idx, bdry in enumerate(self.my_external_nodes):
            self.my_graph.add_edge(infinity_node, bdry, edge_idx)
        if not nx.is_weakly_connected(self.my_graph):
            self.my_graph.remove_node(infinity_node)
            raise ValueError(
                "Every internal vertex should be connected to some boundary vertex")
        self.my_graph.remove_node(infinity_node)

    def __clear_perfect_matching(self):
        """
        clears everything to do with perfect orientation
        """
        self.remove_prop("my_perfect_edge")
        self.my_perfect_matching = None

    def __my_perfect_matching_fix(self, _any_self_loops: bool):
        """
        makes sure the numbers giving my_perfect_edge
        on each vertex induce a perfect matching
        possibly after ignoring some of the boundary vertices
        """
        all_node_names = list(self.my_graph.nodes())
        try:
            special_edge_numbers = {
                z: self.my_graph.nodes[z]["my_perfect_edge"] for z in all_node_names}
        except KeyError:
            self.__clear_perfect_matching()
            return
        self.my_perfect_matching: Optional[Set[Tuple[str, str, int]]] = set()
        already_partnered: Dict[str, str] = {}
        for src_name, tgt_name, edge_key in self.my_graph.edges(keys=True):
            if edge_key == special_edge_numbers[src_name]:
                if src_name in already_partnered:
                    raise ValueError(
                        f"Should not have seen {src_name} as a source yet")
                already_partnered[src_name] = tgt_name
                if self.opposite_colors(src_name,tgt_name) and \
                    already_partnered.get(tgt_name, src_name) != src_name:
                    self.__clear_perfect_matching()
                    return
                if self.same_colors(src_name,tgt_name) and \
                    already_partnered.get(tgt_name, src_name) == src_name:
                    self.__clear_perfect_matching()
                    return
                self.my_perfect_matching.add((src_name, tgt_name, edge_key))
        for src_name in self.my_graph:
            if src_name not in already_partnered and \
                    self.my_graph.nodes[src_name]["is_interior"]:
                self.__clear_perfect_matching()
                break

    def __repair_perfect_matching(self, which_vertices : List[str]) :
        """
        removes the my_perfect_edge from all in which_vertices
        then use the remaining constraint satisfaction problem
        to get new values (now that colors and edges have changed)
        then puts those back in as the new my_perfect_edge
        if no solutions are found, then all perfect matchings are cleared
        """
        if "my_perfect_edge" not in self.my_extra_props:
            self.__clear_perfect_matching()
            return
        which_vertices_2 = (z for z in which_vertices if z in self.my_graph.nodes)
        for cur_node in which_vertices_2:
            try:
                del self.my_graph.nodes[cur_node]["my_perfect_edge"]
            except KeyError:
                pass
        try:
            put_back_iter = self.__remaining_my_perfect_edges()
            put_back = next(put_back_iter)
            found_put_back = True
        except StopIteration:
            found_put_back = False
        except NotImplementedError:
            found_put_back = False
        if found_put_back:
            if set(put_back.keys()) != set(which_vertices):
                self.__clear_perfect_matching()
            else:
                for cur_node in which_vertices :
                    self.__add_props(cur_node,{"my_perfect_edge":put_back[cur_node]})
                self.__my_perfect_matching_fix(False)
        else:
            self.__clear_perfect_matching()

    def __remaining_my_perfect_edges(self) -> Iterator[Dict[str,int]]:
        """
        make the constraint problem for my_perfect_edge
        - each vertex has an integer variable between [0,degree)
        - if an edge connects two opposite colors, then
            r_v==a -> r_w==corresponding(a)
        - if an edge connects two same colors, then
            r_v==a -> r_w!=corresponding(a)
        - for each vertex with a my_perfect_edge prop, that variable
            is pinned down
        """
        for _ in range(0):
            yield {}
        raise NotImplementedError

    def __multi_edge_permutation_add(self, src: str, tgt: str,
                                     src_to_tgt_dict: Optional[Dict[int, int]],
                                     tgt_to_src_dict: Optional[Dict[int, int]],
                                     num_edges: int):
        """
        helper for the initializer/add_nodes with
        the permutation of multi-edge indices from
        the permutation of the half edge numbers
        """
        if num_edges <= 1:
            return
        if src_to_tgt_dict is None and tgt_to_src_dict is None:
            msg = "".join(["When there are multi-edges, a dictionary mapping the edge numbers",
                           f"of the {src} to {tgt} half edges to the {tgt} to {src} half edges"])
            raise ValueError(msg)
        if src_to_tgt_dict is not None and tgt_to_src_dict is not None:
            msg = "".join(["When there are multi-edges,",
                           " exactly one dictionary mapping the edge numbers",
                           f"of the {src} to {tgt} half edges to the {tgt} to {src} half edges",
                           "must be provided. The other direction is the inverse."])
            raise ValueError(msg)

        def make_perm(src_temp: str, tgt_temp: str, which_dict: Dict[int, int]) -> Dict[int, int]:
            """
            check validity of the permutation of the half edge numbers
            which indicates the cyclic orderings
            """
            src_temp_tgt_temp_edge_numbers: List[int] = \
                list(self.my_graph[src_temp][tgt_temp].keys())
            tgt_temp_src_temp_edge_numbers: List[int] = \
                list(self.my_graph[tgt_temp][src_temp].keys())
            if set(which_dict.keys()) != set(src_temp_tgt_temp_edge_numbers):
                msg = " ".join(["When there are multi-edges,",
                               "need a dictionary mapping the edge numbers",
                                f"of the {src} to {tgt} half edges to",
                                "the {tgt} to {src} half edges"])
                raise ValueError(msg)
            if set(which_dict.values()) != tgt_temp_src_temp_edge_numbers:
                msg = " ".join(["When there are multi-edges,",
                               "need a dictionary mapping the edge numbers",
                                f"of the {src} to {tgt} half edges to",
                                "the {tgt} to {src} half edges"])
                raise ValueError(msg)
            return which_dict
        if src_to_tgt_dict is None:
            self.multi_edge_permutation[(tgt, src)] = \
                make_perm(tgt, src, cast(Dict[int, int], tgt_to_src_dict))
        else:
            self.multi_edge_permutation[(src, tgt)] = \
                make_perm(src, tgt, cast(Dict[int, int], src_to_tgt_dict))

    def is_trivalent(self, this_vertex: str) -> bool:
        """
        is this_vertex trivalent
        """
        return len(self.my_graph[this_vertex]) == 3

    def one_color_trivalent(self, all_of_this_color: BiColor) -> bool:
        """
        are all vertices of specified color trivalent
        """
        return all(
            self.get_color(z) != all_of_this_color or self.is_trivalent(z)
            for z in self.my_graph)

    def boundary_connects_to_color(self, all_of_this_color: BiColor) -> bool:
        """
        are boundary vertices connected to only internal vertices of specified color
        no boundary to boundary allowed, no boundary to the opposite color allowed
        """
        for cur_bdry in self.all_bdry_nodes:
            for neighbor in self.my_graph.neighbors(cur_bdry):
                if self.get_color(neighbor) != all_of_this_color:
                    return False
        return True

    def is_bipartite(self) -> bool:
        """
        are all edges connecting internal vertices opposite colors
        """
        for src, tgt in itertools.combinations(self.my_graph, 2):
            if not self.my_graph.nodes[src]["is_interior"] or \
                    not self.my_graph.nodes[tgt]["is_interior"]:
                continue
            if self.nodes_connected(src, tgt) and not self.opposite_colors(src, tgt):
                return False
        return True

    def nodes_connected(self, this_vertex: str, that_vertex: str) -> bool:
        """
        are this_vertex and that_vertex connected by an edge
        """
        num_edges = self.num_connecting_edges(this_vertex, that_vertex)
        return num_edges > 0

    def num_connecting_edges(self, this_vertex: str, that_vertex: str) -> int:
        """
        how many edges connect this_vertex and that_vertex
        """
        try:
            num_edges_connecting = len(self.my_graph[this_vertex][that_vertex])
        except KeyError:
            num_edges_connecting = 0
        return num_edges_connecting

    def __corresponding_multiedge_idx(self, src: str, tgt: str, key_in_src_tgt: int) -> int:
        """
        if we are looking at the half edge labelled key_in_src_tgt'th connecting src to tgt
        then what is the key of the corresponding half edge backwards from tgt to src
        """
        assert 0 <= key_in_src_tgt < self.my_graph.out_degree(src)
        num_edges = self.num_connecting_edges(src, tgt)
        if num_edges == 1:
            all_keys = list(self.my_graph[tgt][src].keys())
            return all_keys[0]
        if num_edges == 0:
            raise ValueError(" ".join([f"There should have been edges from {src} to {tgt}",
                             f"so we could get the one with key {key_in_src_tgt}"]))
        perm = self.multi_edge_permutation.get((src, tgt), None)
        if perm is not None:
            return perm[key_in_src_tgt]
        perm_inv = self.multi_edge_permutation.get((tgt, src), None)
        assert perm_inv is not None
        perm = {}
        to_return = None
        for (key_tgt_src, key_src_tgt) in perm_inv.items():
            if key_src_tgt == key_in_src_tgt:
                to_return = key_tgt_src
            perm[key_src_tgt] = key_tgt_src
        self.multi_edge_permutation[(src, tgt)] = perm
        return cast(int, to_return)

    def get_color(self, this_vertex: str) -> Optional[BiColor]:
        """
        the color of this vertex
        """
        return cast(Optional[BiColor], self.my_graph.nodes()[this_vertex].get("color", None))

    def opposite_colors(self, this_vertex: str, that_vertex: str) -> bool:
        """
        are this_vertex and that_vertex both colored and opposite
        """
        this_color = self.get_color(this_vertex)
        that_color = self.get_color(that_vertex)
        if this_color is None or that_color is None:
            return False
        return this_color != that_color

    def same_colors(self, this_vertex: str, that_vertex: str) -> bool:
        """
        are this_vertex and that_vertex both colored and the same
        """
        this_color = self.get_color(this_vertex)
        that_color = self.get_color(that_vertex)
        if this_color is None or that_color is None:
            return False
        return this_color == that_color

    def _by_edge_number(self, src: str, desired_edge_number: int) -> Tuple[str, int]:
        """
        the target of the half-edge starting at src
        with the desired edge number in the cyclic ordering
        also the corresponding half-edge when looking from target to src
        """
        for tgt in self.my_graph[src]:
            if desired_edge_number in self.my_graph[src][tgt]:
                flipped = self.__corresponding_multiedge_idx(
                    src, tgt, desired_edge_number)
                return tgt, flipped
        raise ValueError(
            f"Nothing with edge number {desired_edge_number} out of {src} found")

    def my_kn_type(self, this_contribution : Optional[str]=None) -> Tuple[int,int]:
        """
        bdry_sinks - bdry_sources + 
            sum_{v in red vertices} (deg(v)-2)
            - sum_{v in green vertices} (deg(v)-2) = 0
        bdry_sinks - bdry_sources = - my_sum
        bdry_sinks + bdry_sources = my_n
        """
        def vertex_contribution(which_vertex : str) -> int:
            """
            contribution of which_vertex to the sum
            """
            cur_color = self.get_color(which_vertex)
            if cur_color is not None:
                cur_degree = cast(int,self.my_graph.out_degree(which_vertex))
                if cur_color == BiColor.RED:
                    my_sum = cur_degree-2
                else:
                    my_sum = -(cur_degree-2)
            else:
                my_sum = 0
            return my_sum
        my_n = len(self.all_bdry_nodes)
        if this_contribution is not None:
            return (vertex_contribution(this_contribution),my_n)
        my_sum = 0
        for vertex_name in self.my_graph.nodes():
            my_sum += vertex_contribution(vertex_name)
        my_sinks = (my_n - my_sum) // 2
        return (my_n-my_sinks,my_n)

    def greedy_shrink(self,
                      name_combiner : Optional[Callable[[str,str],str]] = None,
                      rounds_left : int = 5) -> bool:
        """
        do as many M2 and M3 that remove vertices or edges
        provided no extra data is present
        the name_combiner tells what to do with names on the
        contract edge moves
        """
        if not self.my_extra_props.issubset([]):
            return False
        name_combiner = name_combiner if name_combiner is not None else lambda z1,_: z1
        any_change = False
        all_edges = list(self.my_graph.edges(keys=False))
        for name1,name2 in all_edges:
            if name1 in self.my_graph.nodes() and name2 in self.my_graph.nodes():
                did_remove,_ = self.contract_edge(name1,name2,name_combiner(name1,name2))
            else:
                did_remove = False
            any_change = any_change or did_remove
        all_node_names = list(self.my_graph.nodes())
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
        perform a square move on the quadrilateral with these vertices
        if these names did not meet the pre-conditions for a square move
            provide the explanation as well
        also output if self has changed
        """
        if len(set(four_nodes)) != 4:
            return False, "The 4 need to be distinct"
        all_trivalent = all((self.is_trivalent(cur_node)
                            for cur_node in four_nodes))
        if not all_trivalent:
            return False, "The 4 must be trivalent"
        is_quad = True
        alternating_colors = True

        for cur_node, next_node in zip(four_nodes, four_nodes[1:]):
            is_quad = is_quad and self.nodes_connected(cur_node, next_node)
            alternating_colors = alternating_colors and self.opposite_colors(
                cur_node, next_node)
        is_quad = is_quad and self.nodes_connected(
            four_nodes[3], four_nodes[0])
        alternating_colors = alternating_colors and\
            self.opposite_colors(four_nodes[3], four_nodes[0])
        if not is_quad:
            return False, "The 4 must form a quadrilateral"
        if not alternating_colors:
            return False, "The 4 must alternate colors"

        for cur_node in four_nodes:
            cur_color = self.get_color(cur_node)
            if cur_color is BiColor.RED:
                self.my_graph.nodes()[cur_node]["color"] = BiColor.GREEN
            else:
                self.my_graph.nodes()[cur_node]["color"] = BiColor.RED

        self.__repair_perfect_matching(list(four_nodes))

        return True, "Success"

    # pylint:disable = too-many-return-statements, too-many-statements
    def flip_move(self, this_vertex: str,
                  that_vertex: str,
                  extra_data_transformer:
                  Optional[Callable[[ExtraData, ExtraData, "PlabicGraph"],
                                    Tuple[ExtraData, ExtraData]]] = None)\
            -> Tuple[bool, str]:
        """
        perform a flip move on the edge with these vertices
        if these names did not meet the pre-conditions for a square move
            provide the explanation as well
        also output if self has changed
        """
        if this_vertex == that_vertex:
            return False, "The 2 need to be distinct"
        if not self.is_trivalent(this_vertex) and self.is_trivalent(that_vertex):
            return False, "The 2 must be trivalent"
        if self.num_connecting_edges(this_vertex, that_vertex) != 1:
            msg = "\n".join(["The 2 must have exactly 1 edge connecting them.",
                             "If there were multiple, we don't know which of them",
                             "to collapse then re-expand"])
            return False, msg
        if self.num_connecting_edges(this_vertex, this_vertex) > 0:
            return False, "No self loops for this implementation of flip move to work"
        if self.num_connecting_edges(that_vertex, that_vertex) > 0:
            return False, "No self loops for this implementation of flip move to work"
        first_color = self.get_color(this_vertex)
        second_color = self.get_color(that_vertex)
        if first_color is None or second_color is None or first_color != second_color:
            return False, "The two must be the same color"
        cut_edge_number_this = list(
            self.my_graph[this_vertex][that_vertex].keys())[0]
        cut_edge_number_that = list(
            self.my_graph[that_vertex][this_vertex].keys())[0]

        # if starting from the left hand side picture
        #   if this is the top and that is the bottom
        #       ABCD is upper_left, upper_right, lower_left, lower_right
        #   if this is the bottom and that is the top
        #       ABCD is lower_right, lower_left, upper_right, upper_left
        # if starting from the right hand side picture
        #   if this is the left and that is the right
        #       ABCD is lower_left,upper_left,lower_right,upper_right
        #   if this is the right and that is the left
        #       ABCD is upper_right,lower_right,upper_left,lower_left

        edge_number_a = cut_edge_number_this + 1
        if edge_number_a >= 3:
            edge_number_a -= 3
        vertex_a, halfedge_number_from_a = \
            self._by_edge_number(this_vertex, edge_number_a)

        edge_number_b = cut_edge_number_this - 1
        if edge_number_b < 0:
            edge_number_b += 3
        vertex_b, halfedge_number_from_b = \
            self._by_edge_number(this_vertex, edge_number_b)

        edge_number_d = cut_edge_number_that + 1
        if edge_number_d >= 3:
            edge_number_d -= 3
        vertex_d, halfedge_number_from_d = \
            self._by_edge_number(that_vertex, edge_number_d)

        edge_number_c = cut_edge_number_that - 1
        if edge_number_c < 0:
            edge_number_c += 3
        vertex_c, halfedge_number_from_c = \
            self._by_edge_number(that_vertex, edge_number_c)

        for node_to_go in [this_vertex, that_vertex]:
            adjacencies = list(self.my_graph[node_to_go])
            for adj in adjacencies:
                if self.multi_edge_permutation.get((node_to_go, adj), None) is not None:
                    del self.multi_edge_permutation[(node_to_go, adj)]
                if self.multi_edge_permutation.get((adj, node_to_go), None) is not None:
                    del self.multi_edge_permutation[(adj, node_to_go)]
        if extra_data_transformer is None and \
            self.my_extra_props.issubset(["my_perfect_edge","position"]) and \
            "position" in self.my_extra_props:
            extra_data_transformer = _flip_position_transformer
        if extra_data_transformer is not None:
            old_this_data = cast(ExtraData, self.my_graph.nodes[this_vertex])
            old_that_data = cast(ExtraData, self.my_graph.nodes[that_vertex])
            new_this_data, new_that_data = extra_data_transformer(
                old_this_data, old_that_data, self)
        self.my_graph.remove_node(this_vertex)
        self.my_graph.remove_node(that_vertex)

        self.my_graph.add_node(this_vertex,
                               is_interior=True,
                               color=first_color)
        self.my_graph.add_node(that_vertex,
                               is_interior=True,
                               color=second_color)
        if extra_data_transformer is not None:
            self.__add_props(this_vertex, new_this_data)
            self.__add_props(that_vertex, new_that_data)

        self.my_graph.add_edge(vertex_a, this_vertex, halfedge_number_from_a)
        self.my_graph.add_edge(this_vertex, vertex_a, 2)

        self.my_graph.add_edge(that_vertex, this_vertex, 2)
        self.my_graph.add_edge(this_vertex, that_vertex, 0)

        self.my_graph.add_edge(vertex_c, this_vertex, halfedge_number_from_c)
        self.my_graph.add_edge(this_vertex, vertex_c, 1)

        self.my_graph.add_edge(vertex_b, that_vertex, halfedge_number_from_b)
        self.my_graph.add_edge(that_vertex, vertex_b, 0)

        self.my_graph.add_edge(vertex_d, that_vertex, halfedge_number_from_d)
        self.my_graph.add_edge(that_vertex, vertex_d, 1)

        if vertex_a == vertex_c:
            self.__multi_edge_permutation_add(
                this_vertex, vertex_a,
                {2: halfedge_number_from_a,
                 1: halfedge_number_from_c}, None, 2)
        if vertex_b == vertex_d:
            self.__multi_edge_permutation_add(
                that_vertex, vertex_b,
                {0: halfedge_number_from_b,
                 1: halfedge_number_from_d}, None, 2)

        self.__repair_perfect_matching([this_vertex,that_vertex])

        return True, "Success"

    def remove_bivalent(self, my_bivalent_vertex: str) -> Tuple[bool, str]:
        """
        remove the specified bivalent vertex
        if these names did not meet the pre-conditions of being bivalent
            provide the explanation as well
        also output if self has changed
        """
        should_be_2 = self.my_graph.out_degree(my_bivalent_vertex)
        if should_be_2 != 2:
            return False, f"{my_bivalent_vertex} needs to be bivalent"
        side_1, idx_from_side_1 = self._by_edge_number(my_bivalent_vertex, 0)
        if side_1 == my_bivalent_vertex:
            return False, "No self loops for this implementation of remove bivalent to work"
        side_2, idx_from_side_2 = self._by_edge_number(my_bivalent_vertex, 1)
        if side_2 == my_bivalent_vertex:
            return False, "No self loops for this implementation of remove bivalent to work"
        if self.get_color(side_1) is None and self.get_color(side_2) is None:
            return False, "Cannot cause two boundary nodes to connect directly"
        num_side_1_to_side_2 = self.num_connecting_edges(side_1, side_2)
        if num_side_1_to_side_2 > 1:
            old_dict_12 = self.multi_edge_permutation.get(
                (side_1, side_2), None)
            old_dict_21 = self.multi_edge_permutation.get(
                (side_2, side_1), None)
        if num_side_1_to_side_2 == 1:
            all_keys_12 = list(self.my_graph[side_1][side_2].keys())
            all_keys_21 = list(self.my_graph[side_2][side_1].keys())
            old_dict_12 = {all_keys_12[0]: all_keys_21[0]}
            old_dict_21 = None
        if num_side_1_to_side_2 == 0:
            old_dict_12 = None
            old_dict_21 = None
        for adj in [side_1, side_2]:
            if self.multi_edge_permutation.get((my_bivalent_vertex, adj), None) is not None:
                del self.multi_edge_permutation[(my_bivalent_vertex, adj)]
            if self.multi_edge_permutation.get((adj, my_bivalent_vertex), None) is not None:
                del self.multi_edge_permutation[(adj, my_bivalent_vertex)]
        self.my_graph.remove_node(my_bivalent_vertex)
        self.my_graph.add_edge(side_1, side_2, idx_from_side_1)
        self.my_graph.add_edge(side_2, side_1, idx_from_side_2)
        #pylint:disable=possibly-used-before-assignment
        if old_dict_12 is not None:
            old_dict_12[idx_from_side_1] = idx_from_side_2
            self.multi_edge_permutation[(side_1, side_2)] = old_dict_12
        #pylint:disable=possibly-used-before-assignment
        if old_dict_21 is not None:
            old_dict_21[idx_from_side_2] = idx_from_side_1
            self.multi_edge_permutation[(side_2, side_1)] = old_dict_21

        self.__repair_perfect_matching([side_1,side_2])

        return True, "Success"

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
        num_edges_connecting = self.num_connecting_edges(
            this_vertex, that_vertex)
        if num_edges_connecting != 1:
            return False, "There needs to be exactly one edge to split specified by these vertices"
        if desired_name in self.my_graph.nodes():
            return False, "The new vertex name should be distinct from those present"
        key_this_that_all: List[int] = list(
            self.my_graph[this_vertex][that_vertex].keys())
        key_this_that = key_this_that_all[0]
        key_that_this = self.__corresponding_multiedge_idx(
            this_vertex, that_vertex, key_this_that)

        if extra_data_transformer is None and \
            self.my_extra_props.issubset(["my_perfect_edge","position"]) and \
            "position" in self.my_extra_props:
            extra_data_transformer = _insert_bivalent_position_transformer
        if extra_data_transformer is not None:
            old_this_data = cast(ExtraData, self.my_graph.nodes[this_vertex])
            old_that_data = cast(ExtraData, self.my_graph.nodes[that_vertex])
            new_props = extra_data_transformer(
                old_this_data, old_that_data, self)
        self.my_graph.add_node(
            desired_name, is_interior=True, color=desired_color)
        if extra_data_transformer is not None:
            self.__add_props(desired_name, new_props)
        self.my_graph.remove_edge(this_vertex, that_vertex)
        self.my_graph.add_edge(this_vertex, desired_name, key_this_that)
        self.my_graph.add_edge(desired_name, this_vertex, 0)
        self.my_graph.add_edge(that_vertex, desired_name, key_that_this)
        self.my_graph.add_edge(desired_name, that_vertex, 1)

        self.__repair_perfect_matching([this_vertex,that_vertex,desired_name])

        return True, "Success"

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
        if this_vertex == that_vertex:
            return False, "No self loops for this implementation of contract edge to work"
        this_color = self.get_color(this_vertex)
        that_color = self.get_color(that_vertex)
        if this_color is None or that_color is None:
            return False, "The two vertices for the collapsing edge must be internal"
        if this_color != that_color:
            return False, "They must be the same color, to collapse into 1"
        num_edges_connecting = self.num_connecting_edges(
            this_vertex, that_vertex)
        if num_edges_connecting != 1:
            return False,\
                "There needs to be exactly one edge to collapse specified by these vertices"
        if combined_name in self.my_graph.nodes() and \
                combined_name not in [this_vertex, that_vertex]:
            return False, "The new vertex name should be distinct from those present"
        all_adjacencies: List[str] = list(self.my_graph[this_vertex])
        all_adjacencies.extend(self.my_graph[that_vertex])
        connects_to_this = [False] * len(all_adjacencies)
        connects_to_that = [False] * len(all_adjacencies)
        will_multiconnect = []
        for idx, cur_adj in enumerate(all_adjacencies):
            num_to_this = self.num_connecting_edges(cur_adj, this_vertex)
            connects_to_this[idx] = num_to_this > 0
            num_to_that = self.num_connecting_edges(cur_adj, that_vertex)
            connects_to_that[idx] = num_to_that > 0
            if num_to_this+num_to_that > 1:
                will_multiconnect.append(cur_adj)
        key_this_that_all: List[int] = list(
            self.my_graph[this_vertex][that_vertex].keys())
        key_this_that = key_this_that_all[0]
        out_this = self.my_graph.out_degree(this_vertex)
        key_that_this = self.__corresponding_multiedge_idx(
            this_vertex, that_vertex, key_this_that)
        out_that = self.my_graph.out_degree(that_vertex)
        order_around_combined: List[Tuple[bool, int, Tuple[str, int]]] = []
        order_around_combined.extend(((True, z, self._by_edge_number(this_vertex, z))
                                      for z in range(key_this_that+1, out_this)))
        order_around_combined.extend(((True, z, self._by_edge_number(this_vertex, z))
                                      for z in range(0, key_this_that)))
        order_around_combined.extend(((False, z, self._by_edge_number(that_vertex, z))
                                      for z in range(key_that_this+1, out_that)))
        order_around_combined.extend(((False, z, self._by_edge_number(that_vertex, z))
                                      for z in range(0, key_that_this)))

        for node_to_go in [this_vertex, that_vertex]:
            adjacencies = list(self.my_graph[node_to_go])
            for adj in adjacencies:
                if self.multi_edge_permutation.get((node_to_go, adj), None) is not None:
                    del self.multi_edge_permutation[(node_to_go, adj)]
                if self.multi_edge_permutation.get((adj, node_to_go), None) is not None:
                    del self.multi_edge_permutation[(adj, node_to_go)]
        if extra_data_transformer is None and \
            self.my_extra_props.issubset(["my_perfect_edge","position"]) and \
            "position" in self.my_extra_props:
            extra_data_transformer = _contract_edge_position_transformer
        if extra_data_transformer is not None:
            old_this_data = cast(ExtraData, self.my_graph.nodes[this_vertex])
            old_that_data = cast(ExtraData, self.my_graph.nodes[that_vertex])
            new_props = extra_data_transformer(
                old_this_data, old_that_data, self)
        self.my_graph.remove_node(this_vertex)
        self.my_graph.remove_node(that_vertex)

        self.my_graph.add_node(
            combined_name, is_interior=True, color=this_color)
        if extra_data_transformer is not None:
            self.__add_props(combined_name, new_props)
        for key_from_combined, (_origin_this, _origin_key,
                                (tgt, half_edge_from_tgt)) in\
                enumerate(order_around_combined):
            self.my_graph.add_edge(combined_name, tgt, key_from_combined)
            self.my_graph.add_edge(tgt, combined_name,
                                   half_edge_from_tgt)
            if tgt in will_multiconnect:
                if (combined_name, tgt) in self.multi_edge_permutation:
                    self.multi_edge_permutation[(combined_name, tgt)][key_from_combined] =\
                        half_edge_from_tgt
                else:
                    self.multi_edge_permutation[(combined_name, tgt)] = {}
                    self.multi_edge_permutation[(combined_name, tgt)][key_from_combined] =\
                        half_edge_from_tgt

        self.__repair_perfect_matching([combined_name])

        return True, "Success"

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
        this_color = self.get_color(this_vertex)
        if this_color is None:
            return False, "The vertex for splitting must be internal"
        if split_name_1 == split_name_2:
            return False, "The new vertex names should be distinct from each other"
        if split_name_1 in self.my_graph.nodes() or \
                split_name_2 in self.my_graph.nodes():
            if split_name_1 == this_vertex and split_name_2 not in self.my_graph.nodes():
                pass
            elif split_name_2 == this_vertex and split_name_1 not in self.my_graph.nodes():
                pass
            else:
                return False, "The new vertex names should be distinct from those present"
        this_out_degree = self.my_graph.out_degree(this_vertex)
        if split_bounds[0] <= split_bounds[1]:
            go_to_1 = list(range(split_bounds[0], split_bounds[1]+1))
            go_to_2 = list(range(split_bounds[1]+1, this_out_degree)) +\
                list(range(0, split_bounds[0]))
        else:
            go_to_2 = list(range(split_bounds[1], split_bounds[0]+1))
            go_to_1 = list(range(split_bounds[0]+1, this_out_degree)) +\
                list(range(0, split_bounds[1]))

        go_to_1_full : List[Tuple[int,str,int]] = []
        seen_on_1 : Set[str] = set()
        will_multiconnect_1 : Set[str] = set()
        for on_1_key, on_this_key in enumerate(go_to_1):
            tgt, on_tgt_key = self._by_edge_number(this_vertex, on_this_key)
            go_to_1_full.append((on_1_key, tgt, on_tgt_key))
            if tgt in seen_on_1:
                will_multiconnect_1.add(tgt)
            seen_on_1.add(tgt)
        go_to_2_full : List[Tuple[int,str,int]] = []
        seen_on_2 : Set[str] = set()
        will_multiconnect_2 : Set[str] = set()
        for on_2_key, on_this_key in enumerate(go_to_2):
            tgt, on_tgt_key = self._by_edge_number(this_vertex, on_this_key)
            go_to_2_full.append((on_2_key, tgt, on_tgt_key))
            if tgt in seen_on_2:
                will_multiconnect_2.add(tgt)
            seen_on_2.add(tgt)

        adjacencies = list(self.my_graph[this_vertex])
        for adj in adjacencies:
            if self.multi_edge_permutation.get((this_vertex, adj), None) is not None:
                del self.multi_edge_permutation[(this_vertex, adj)]
            if self.multi_edge_permutation.get((adj, this_vertex), None) is not None:
                del self.multi_edge_permutation[(adj, this_vertex)]
        if extra_data_transformer is None and \
            self.my_extra_props.issubset(["my_perfect_edge","position"]) and \
            "position" in self.my_extra_props:
            extra_data_transformer = functools.partial(
                _split_vertex_position_transformer,this_vertex,(seen_on_1,seen_on_2))
        if extra_data_transformer is not None:
            old_this_data = cast(ExtraData, self.my_graph.nodes[this_vertex])
            new_props_1, new_props_2 = extra_data_transformer(
                old_this_data, self)
        self.my_graph.remove_node(this_vertex)

        self.my_graph.add_node(
            split_name_1, is_interior=True, color=this_color)
        self.my_graph.add_node(
            split_name_2, is_interior=True, color=this_color)
        if extra_data_transformer is not None:
            self.__add_props(split_name_1, new_props_1)
            self.__add_props(split_name_2, new_props_2)
        out_on_1 = len(go_to_1_full)
        out_on_2 = len(go_to_2_full)
        self.my_graph.add_edge(split_name_1, split_name_2, out_on_1)
        self.my_graph.add_edge(split_name_2, split_name_1, out_on_2)
        for on_1_key, tgt, on_tgt_key in go_to_1_full:
            self.my_graph.add_edge(split_name_1, tgt, on_1_key)
            self.my_graph.add_edge(tgt, split_name_1, on_tgt_key)
            if tgt in will_multiconnect_1:
                if (split_name_1, tgt) not in self.multi_edge_permutation:
                    self.multi_edge_permutation[(split_name_1, tgt)] = {}
                self.multi_edge_permutation[(
                    split_name_1, tgt)][on_1_key] = on_tgt_key
        for on_2_key, tgt, on_tgt_key in go_to_2_full:
            self.my_graph.add_edge(split_name_2, tgt, on_2_key)
            self.my_graph.add_edge(tgt, split_name_2, on_tgt_key)
            if tgt in will_multiconnect_2:
                if (split_name_2, tgt) not in self.multi_edge_permutation:
                    self.multi_edge_permutation[(split_name_2, tgt)] = {}
                self.multi_edge_permutation[(
                    split_name_2, tgt)][on_2_key] = on_tgt_key

        self.__repair_perfect_matching([split_name_1,split_name_2])

        return True, "Success"

    def __operad_compose_helper(self,which_internal_disk : int,
                                   other:PlabicGraph,
                                   glued_vertices : List[str],
                                   force_it : bool = False)\
                                    -> Tuple[bool,Callable[[Point],Point]]:
        """
        find a transformation
        taking the positions of the glued vertices in self
        to the corresponding positions in other
        does that transformation to self.circles_config
        """
        try:
            all_positions_self = {
                    z: self.my_graph.nodes[z]["position"] for z in glued_vertices}
            all_positions_other = {
                    z: other.my_graph.nodes[z]["position"] for z in glued_vertices}
        except KeyError:
            # something didn't have position so can't trust circles_config
            # or positions in the operad composition
            # circles_config is removed here
            # the False output tells operad_compose to clear position property
            self.circles_config = None
            return False,lambda z: z
        if self.circles_config is not None and other.circles_config is not None:
            supposed_transformation = \
                self.circles_config.do_transformation(which_internal_disk,other.circles_config)
            transformed_self = {z1:supposed_transformation(cast(Point,z2))
                                for (z1,z2) in all_positions_self.items()}
            if transformed_self == all_positions_other:
                # the transformations of the framed circles accomplished the task
                # the circles_config has already changed
                # the points will all change later in operad_compose
                return True,supposed_transformation
            # the transformations of the framed circles failed the task
            # circles_config is removed here
            # the False output tells operad_compose to clear position property
            self.circles_config = None
            return False,lambda z: z
        if all_positions_self == all_positions_other:
            # the points lined up but we didn't have both circles_config for self and other
            # so we can't have one in the composition
            # but we still have the positions of the vertices in the operad composition
            self.circles_config = None
            return True,lambda z: z
        if not force_it:
            # the points did not line up and we didn't have both circles_config for self and other
            # so we can't have one in the composition
            # and we aren't forcing to find an affine transformation that does the task without
            # the help of FramedDisKConfig
            self.circles_config = None
            return False,lambda z: z
        raise NotImplementedError

    def radius_unoccupied(self, center : Point, ignored_vertices : List[str]) -> Tuple[bool,float]:
        """
        the ball of the returned radius around center
        is only allowed to have ignored_vertices within it
        """
        def point_like(maybe_point) -> bool:
            """
            is p in R^2
            """
            return isinstance(maybe_point,tuple) and len(maybe_point)==2 and\
                isinstance(maybe_point[0],numbers.Real) and isinstance(maybe_point[1],numbers.Real)
        all_are_points = True
        for node_name in self.my_graph.nodes():
            old_position = self.my_graph.nodes[node_name].get("position",None)
            if not point_like(old_position):
                all_are_points = False
                break
        if not all_are_points:
            return False,0.0
        min_radius_sq : Optional[float] = None
        for node_name in self.my_graph.nodes():
            if node_name in ignored_vertices:
                continue
            cur_position = cast(Point,self.my_graph.nodes[node_name]["position"])
            cur_radius_sq = (cur_position[0]-center[0])**2 + (cur_position[1]-center[1])**2
            if min_radius_sq is None or cur_radius_sq<min_radius_sq:
                min_radius_sq = cur_radius_sq
        if min_radius_sq is None:
            if self.circles_config is None:
                return False, 1.0
            (disk_center, disk_radius, _) = self.circles_config.outer_circle
            center_to_disk_center = (disk_center[0]-center[0])**2 + (disk_center[1]-center[1])**2
            if center_to_disk_center<disk_radius**2:
                return True, disk_radius-sqrt(center_to_disk_center)
            return False, 0.0
        return True,sqrt(min_radius_sq)

    def coordinate_transform(self,transform : Callable[[Point],Point],
                             circle_config_invalidates : bool = True) -> bool:
        """
        does transform on all positions of all nodes
        """
        def point_like(maybe_point) -> bool:
            """
            is p in R^2
            """
            return isinstance(maybe_point,tuple) and len(maybe_point)==2 and\
                isinstance(maybe_point[0],numbers.Real) and isinstance(maybe_point[1],numbers.Real)
        all_are_points = True
        for node_name in self.my_graph.nodes():
            old_position = self.my_graph.nodes[node_name].get("position",None)
            if not point_like(old_position):
                # print(f"{node_name} gave {old_position}")
                all_are_points = False
                break
        if all_are_points:
            if circle_config_invalidates:
                self.circles_config = None
            for node_name in self.my_graph.nodes():
                old_position = self.my_graph.nodes[node_name]["position"]
                new_position = transform(old_position)
                self.my_graph.nodes[node_name]["position"] = new_position
        return all_are_points

    def operad_compose(self, other: PlabicGraph, which_internal_disk: int) -> Tuple[bool, str]:
        """
        substitute other in on the i'th internal circle of self
            that is self has at least that many in self.my_internal_bdry
            and self.my_internal_bdry[which_internal_disk] are matched up
            with the other.my_external_nodes
        if anything did not meet the pre-conditions
            provide the explanation as well
        also output if self has changed
        """
        if which_internal_disk < 0:
            return False, "Which internal disk is a natural number"
        if self.num_interior_circles <= which_internal_disk:
            return False, "Not enough internal disks"
        relevant_internal = self.my_internal_bdry[which_internal_disk]

        if not cyclic_equal(relevant_internal, other.my_external_nodes):
            return False, "The boundary vertices that should be glued did not match up"

        nodes_overlapping = set(self.my_graph.nodes())
        nodes_overlapping = nodes_overlapping.intersection(
            set(other.my_graph.nodes()))
        if nodes_overlapping != set(relevant_internal):
            return False, "The only overlap of vertex names should be the ones that are glued"

        for replaced_node in relevant_internal:
            self_tgt, _ = self._by_edge_number(replaced_node, 0)
            # pylint:disable = protected-access
            other_tgt, _ = other._by_edge_number(replaced_node, 0)
            if self.get_color(self_tgt) is None or \
                    other.get_color(other_tgt) is None:
                return False,\
                    "Boundary vertices should connect to internal vertices not directly to boundary"
        if "position" in self.my_extra_props and "position" in other.my_extra_props:
            found_transformation, my_transform = \
                self.__operad_compose_helper(which_internal_disk, other,relevant_internal)
            clear_position = not found_transformation
            if found_transformation:
                self.coordinate_transform(my_transform,False)
        else:
            clear_position = True

        self.all_bdry_nodes = self.my_external_nodes
        new_internal_bdry: List[List[str]] = \
            [[]]*(self.num_interior_circles-1+other.num_interior_circles)
        for idx_l in range(which_internal_disk):
            new_internal_bdry[idx_l] = self.my_internal_bdry[idx_l]
            self.all_bdry_nodes.extend(new_internal_bdry[idx_l])
        idx_l = which_internal_disk
        for idx_r in range(other.num_interior_circles):
            new_internal_bdry[idx_l+idx_r] = other.my_internal_bdry[idx_r]
            self.all_bdry_nodes.extend(new_internal_bdry[idx_l+idx_r])
        idx_l = other.num_interior_circles
        for idx_r in range(which_internal_disk, self.num_interior_circles):
            new_internal_bdry[idx_l+idx_r] = self.my_internal_bdry[idx_r]
            self.all_bdry_nodes.extend(new_internal_bdry[idx_l+idx_r])
        self.my_internal_bdry = new_internal_bdry
        self.num_interior_circles += other.num_interior_circles - 1

        for replaced_node in relevant_internal:
            tgt, _ = self._by_edge_number(replaced_node, 0)
            self.my_graph.add_edge(replaced_node, tgt, 1)
            self.my_graph.remove_edge(replaced_node, tgt, 0)
            self.my_graph.nodes[replaced_node]["is_interior"] = True
            self.my_graph.nodes[replaced_node]["color"] = BiColor.RED
        self.my_graph = nx.compose(other.my_graph, self.my_graph)

        self.multi_edge_permutation.update(other.multi_edge_permutation)
        for replaced_node in relevant_internal:
            bivalent_removed, msg = self.remove_bivalent(replaced_node)
            if not bivalent_removed:
                new_msg = " ".join([f"Got the message {msg} when trying",
                                   f"to remove {replaced_node} as a bivalent vertex.",
                                    "It should have been bivalent because in the constituents",
                                    "it was univalent to distinct targets in the two pieces",
                                    "then the composition would make it bivalent"])
                raise ValueError(new_msg)

        if clear_position:
            self.remove_prop("position")
        # clears everything to do with perfect orientation
        # could fix the perfect orientation instead
        # along this move
        self.__clear_perfect_matching()

        return True, "Success"

    def bdry_to_bdry(self, starting_bdry_node: str,
                     care_about_lollipop : bool = True) -> Tuple[Optional[BiColor], str]:
        """
        starting at a boundary vertex follow the rules of the road until another
        boundary vertex (possibly the same one)
        the color gives the decorated permutation
        because fixed points are colored by the color of the lollipop
        """
        was_at_this_before = None
        via_edge_keyed = None
        at_this_node = starting_bdry_node
        turn_around_colors: List[BiColor] = []
        path_length = 0
        while True:
            next_node = self.follow_rules_of_road(
                at_this_node, was_at_this_before, via_edge_keyed)
            if next_node is None:
                if at_this_node == starting_bdry_node:
                    if len(turn_around_colors) == 0 or not care_about_lollipop:
                        return None, at_this_node
                    if path_length == 2:
                        return turn_around_colors[0], at_this_node
                    msg = " ".join(["There was a collapsible tree",
                                    "but couldn't tell what decoration",
                                    "it collapsed to"])
                    raise NotImplementedError(msg)
                return None, at_this_node
            path_length += 1
            if was_at_this_before is not None and next_node[1] == was_at_this_before:
                turn_around_colors.append(self.get_color(at_this_node))
            was_at_this_before = at_this_node
            via_edge_keyed, at_this_node = next_node

    def follow_rules_of_road(self, at_this_node: str,
                             was_at_this_before: Optional[str] = None,
                             via_edge_keyed: Optional[int] = None) -> Optional[Tuple[int, str]]:
        """
        if we are at some vertex and came from some other node, where should we go next
        if we are following the rule that we turn maximally right/left
        depending on the color of where we are and stop if we are on a boundary vertex
        """
        if was_at_this_before is None:
            if not at_this_node in self.all_bdry_nodes:
                raise ValueError(
                    "If don't have a previous vertex, we must be at a boundary vertex")
            return_name, _ = self._by_edge_number(at_this_node, 0)
            return 0, return_name
        if via_edge_keyed is None:
            only_one = self.num_connecting_edges(
                was_at_this_before, at_this_node) == 1
            if only_one:
                key_this_that_all: List[int] =\
                    list(self.my_graph[was_at_this_before]
                         [at_this_node].keys())
                via_edge_keyed = key_this_that_all[0]
            else:
                msg = " ".join(["If there was a previous vertex,",
                                "then there had to have been a",
                                "previous key for the half edge going from",
                                "previous to current"])
                raise ValueError(msg)
        my_color = self.get_color(at_this_node)
        if my_color is None:
            return None
        prev_edge_number = \
            self.__corresponding_multiedge_idx(
                was_at_this_before, at_this_node, via_edge_keyed)
        num_incident_edges = self.my_graph.out_degree(at_this_node)
        if my_color == BiColor.GREEN:
            next_edge_number = prev_edge_number + 1
            if next_edge_number >= num_incident_edges:
                next_edge_number -= num_incident_edges
        else:
            next_edge_number = prev_edge_number - 1
            if next_edge_number < 0:
                next_edge_number += num_incident_edges
        for tgt in self.my_graph[at_this_node]:
            if next_edge_number in self.my_graph[at_this_node][tgt]:
                return next_edge_number, tgt
        raise ValueError(
            f"There should be an edge numbered {next_edge_number} on {at_this_node}")

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
        color_dict = {BiColor.RED: red_nodes, BiColor.GREEN: green_nodes}

        def name_to_color(name: str) -> str:
            """
            convert from node name to the color to draw it as
            """
            return color_dict.get(
                cast(BiColor, self.get_color(name)), bdry_nodes)

        all_node_names = list(self.my_graph.nodes())
        all_colors = [name_to_color(z) for z in all_node_names]
        if "position" in self.my_extra_props:
            try:
                all_positions = {
                    z: self.my_graph.nodes[z]["position"] for z in all_node_names}
            except KeyError:
                self.remove_prop("position")
                all_positions = None
        else:
            all_positions = None

        something_transparent = "#0f0f0f00"
        if overridden_arrow_orientation is None:
            draw_arrowheads = False
            try:
                special_edge_numbers = {
                    z: self.my_graph.nodes[z]["my_perfect_edge"] for z in all_node_names}
                edge_is_special: List[Tuple[str, str, bool]] = \
                    [(u, v, k == special_edge_numbers[u])
                    for u, v, k in self.my_graph.edges(keys=True)]
                if draw_oriented_if_perfect or not self.is_bipartite():
                    draw_arrowheads = True

                    def keep_this_arrow(src_name: str, tgt_name: str, is_special: bool) -> bool:
                        """
                        whether this arrow is kept or not
                        """
                        src_color = self.get_color(src_name)
                        if src_color is None:
                            tgt_color = self.get_color(tgt_name)
                            if tgt_color is None:
                                # boundary to boundary edge
                                # orientation not determined by color
                                # of endpoints
                                return (src_name < tgt_name) ^ is_special
                            return (tgt_color == BiColor.GREEN) ^ (not is_special)
                        return (src_color == BiColor.RED) ^ (not is_special)
                    edge_colors = [oriented_arrows if keep_this_arrow(u, v, is_special)
                                else something_transparent for u, v, is_special in edge_is_special]
                else:
                    edge_colors = [unoriented_arrows_perfect if is_special
                                else unoriented_arrows_imperfect
                                for u, v, is_special in edge_is_special]
            except KeyError:
                edge_colors = None
        else:
            edge_colors = [oriented_arrows if overridden_arrow_orientation(u, v, k)
                        else something_transparent
                        for u, v, k in self.my_graph.edges(keys=True)]
            draw_arrowheads = None
        if draw_arrowheads is not None:
            nx.draw(self.my_graph, pos=all_positions,
                    node_color=all_colors, edge_color=edge_colors,
                    arrows=draw_arrowheads, with_labels=show_node_names)
        else:
            nx.draw(self.my_graph, pos=all_positions,
                    node_color=all_colors, edge_color=edge_colors,
                    with_labels=show_node_names)
        if self.circles_config is not None:
            self.circles_config.draw(on_top_of_existing=True,
                show_as_well=False,
                draw_framing=False,
                external_circle_color=external_circle_color,
                internal_circles_color=internal_circles_color)
        plt.draw()
        if show_as_well:
            plt.show()

def _flip_position_transformer(data_1 : ExtraData,
                                data_2 : ExtraData,
                                surrounding_plabic : PlabicGraph) -> Tuple[ExtraData, ExtraData]:
    """
    the default extra data transformer for flip move
    """
    d1_pos : Point = data_1["position"]
    d2_pos : Point = data_2["position"]
    halfway_btw = ((d1_pos[0]+d2_pos[0])/2.0,(d1_pos[1]+d2_pos[1])/2.0)
    #pylint:disable=pointless-string-statement
    """make the segment perpendicular to the line d1-d2
    passing through the midpoint
    the new positions are the endpoints of that segment
    use the rest of the positions in _p to determine how long the segment is
    the same length as the original segment unless there are other points in the way"""
    found_rad, rad_unoccupied = surrounding_plabic.radius_unoccupied(halfway_btw, [])
    displacement = (d1_pos[0]-d2_pos[0],d1_pos[1]-d2_pos[1])
    dist_d1_d2 = sqrt(displacement[0]**2+displacement[1]**2)
    _radius_to_use = min(dist_d1_d2/2.0, rad_unoccupied) if found_rad else dist_d1_d2/2.0
    raise NotImplementedError("__flip_position_transformer")

def _insert_bivalent_position_transformer(data_1 : ExtraData,
                                           data_2 : ExtraData,
                                           _surrounding_plabic : PlabicGraph) -> ExtraData:
    """
    the default extra data transformer for insert bivalent move
    """
    d1_pos = data_1["position"]
    d2_pos = data_2["position"]
    halfway_btw = ((d1_pos[0]+d2_pos[0])/2.0,(d1_pos[1]+d2_pos[1])/2.0)
    return {"position":halfway_btw}

def _contract_edge_position_transformer(data_1 : ExtraData,
                                         data_2 : ExtraData,
                                         _surrounding_plabic : PlabicGraph) -> ExtraData:
    """
    the default extra data transformer for contract edge move
    """
    d1_pos = data_1["position"]
    d2_pos = data_2["position"]
    halfway_btw = ((d1_pos[0]+d2_pos[0])/2.0,(d1_pos[1]+d2_pos[1])/2.0)
    return {"position":halfway_btw}

def _split_vertex_position_transformer(my_name : str,_nhbrs_on_each : Tuple[Set[str]],
                                       data_1 : ExtraData,
                                        surrounding_plabic : PlabicGraph)\
                                            -> Tuple[ExtraData, ExtraData]:
    """
    the default extra data transformer for split vertex move
    """
    # find the radius around d1_pos that is not occupied by any other vertices
    # the 2 vertices will need to be somewhere within that ball
    d1_pos = data_1["position"]
    found_rad, rad_unoccupied = surrounding_plabic.radius_unoccupied(d1_pos, [my_name])
    if not found_rad:
        raise NotImplementedError("__split_vertex_position_transformer")
    _position_1 = (d1_pos[0]-rad_unoccupied*1.0/3.0,d1_pos[1])
    _position_2 = (d1_pos[0]+rad_unoccupied*1.0/3.0,d1_pos[1])
    raise NotImplementedError("__split_vertex_position_transformer")
    # use the nhbrs_on_each and their positions to determine which position
    #   for which new vertex
    # return ({"position":position_1},{"position":position_2})
