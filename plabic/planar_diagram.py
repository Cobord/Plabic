"""
diagrams used in the sense of cluster algebra
"""

from enum import Enum, auto
from typing import Iterator, List, Set, Tuple, TypeVar, Union, Callable, Optional
from typing import runtime_checkable, Protocol, Any, cast
import itertools
from functools import reduce
import networkx as nx
import matplotlib.pyplot as plt

T = TypeVar("T")

class Arithmetic(Protocol):
    """
    can be multiplied and added
    """
    def __mul__(self : T,other : T) -> T: ...
    def __add__(self : T,other : T) -> T: ...

NT = TypeVar("NT",bound=Arithmetic)

WeightType = Union[float, NT]
EdgeData = Tuple[int, int, WeightType]
Nat = int

class ChipType(Enum):
    """
    which direction the line on a chip goes
    """
    UP = auto()
    DOWN = auto()
    FLAT = auto()


ChipWord = List[Tuple[ChipType, Nat]]

# pylint:disable=line-too-long,invalid-name,too-many-locals,no-member,too-many-branches,too-few-public-methods


class PlanarNetwork:
    """
    built on top of a directed graph
    """

    @staticmethod
    def from_chip_word(chip_word: ChipWord, weights: List[WeightType]) -> List[List[EdgeData]]:
        """
        edge data as in initializer
        the initializer doesn't take this because we also want to
        present the diagram as it would be drawn
        with multiple elementary pieces occuring in the same stage
        provided they concern different indices
        """
        if len(weights) != len(chip_word):
            raise ValueError("The number of weights and chips must match")
        ret_val = []
        for ((ct, idx), wt) in zip(chip_word, weights):
            if ct == ChipType.UP:
                cur_edge = (idx, idx+1, wt)
            elif ct == ChipType.DOWN:
                cur_edge = (idx-1, idx, wt)
            else:
                cur_edge = (idx, idx, wt)
            ret_val.append([cur_edge])
        return ret_val

    #pylint:disable=too-many-arguments
    def __init__(self, n: Nat, *, edge_list: Optional[List[List[EdgeData]]] = None,
                 chip_word: Optional[ChipWord] = None, chip_weights: Optional[List[WeightType]] = None,
                 multiplicative_identity: WeightType = 1.0,
                 additive_identity: WeightType = 0.0,
                 totally_connected: bool = True):
        """
        a planar diagram with n horizontal lines
        """
        self._underlying_graph = nx.DiGraph()
        self.multiplicative_identity = multiplicative_identity
        self.additive_identity = additive_identity
        self.totally_connected = totally_connected
        for x in range(1, n+1):
            self._underlying_graph.add_node(f"{x}_0",position=(0,x))
        sink_number = {x: 0 for x in range(1, n+1)}
        if edge_list is None:
            if chip_word is not None and chip_weights is not None:
                edge_list = PlanarNetwork.from_chip_word(
                    chip_word, chip_weights)
            else:
                raise ValueError(
                    "Either edge_list must be given or both chip_word and chip_weights")
        self.chip_type: Optional[ChipWord] = [] if chip_word is None else chip_word
        # edge_list.append([(idx,idx,multiplicative_identity) for idx in range(1,n+1)])
        for stage_num,stage in enumerate(edge_list):
            already_incremented = []
            for edge in stage:
                from_node, to_node, cur_weight = edge
                if chip_word is None and self.chip_type is not None:
                    if from_node-to_node == 1:
                        self.chip_type.append((ChipType.DOWN, to_node))
                    elif to_node - from_node == 1:
                        self.chip_type.append((ChipType.UP, to_node))
                    elif to_node == from_node:
                        self.chip_type.append((ChipType.FLAT, to_node))
                    else:
                        self.chip_type = None
                if from_node in already_incremented:
                    from_number = sink_number[from_node]-1
                else:
                    from_number = sink_number[from_node]
                if to_node in already_incremented:
                    to_number = sink_number[to_node]
                else:
                    to_number = sink_number[to_node]+1
                    sink_number[to_node] += 1
                    self._underlying_graph.add_edge(f"{to_node}_{to_number-1}",
                                                    f"{to_node}_{to_number}", weight=multiplicative_identity)
                    already_incremented.append(to_node)
                to_position = (stage_num+1,to_node)
                self._underlying_graph.add_node(f"{to_node}_{to_number}",position=to_position)
                self._underlying_graph.add_edge(f"{from_node}_{from_number}",
                                                f"{to_node}_{to_number}", weight=cur_weight)
        self.sink_number = sink_number

    def path_weight(self, path_nodes: List[Any]) -> WeightType:
        """
        multiply the weights for the path that connects all the path_nodes
        """
        path_edge_weights : Iterator[WeightType] = (self._underlying_graph[path_nodes[i]][path_nodes[i+1]]['weight']
                      for i in range(len(path_nodes)-1))
        return reduce(lambda acc, y: acc*y, (z for z in path_edge_weights), self.multiplicative_identity)

    def weight_matrix(self, my_i: Nat, my_j: Nat) -> WeightType:
        """
        the sum of the path weights over all simple paths connecting
        the source node on horizontal line i
        the sink node on horizontal line j
        """
        source_node = f"{my_i}_0"
        sink_node = f"{my_j}_{self.sink_number[my_j]}"
        weight_ij = self.additive_identity
        for path_nodes in nx.all_simple_paths(self._underlying_graph, source_node, sink_node):
            cur_contrib = self.path_weight(path_nodes)
            weight_ij += cur_contrib
        return weight_ij

    def vertex_disjoint_collection(self, i_set: Set[Nat], j_set: Set[Nat]) -> Iterator[WeightType]:
        """
        iterate through all systems of paths
        that connect the sources with i_set to the sinks with j_set
        suppose there are k of each so we get k paths each time
        skip if they are not vertex disjoint
        but if they are vertex disjoint give the product of the path weights
        """
        # start off by finding the possibilities of how i_set and j_set match up
        # so that (i,j) in collection means that particular i in i_set connects to j in j_set
        if len(i_set) == 0 or len(j_set) == 0:
            raise ValueError("The sets of sources and sinks must be nonempty")
        relevant_pairings = cast(List[Set[Tuple[Nat,Nat]]],[pairings_nc(i_set, j_set)])
        # relevant_pairings = pairings(i_set,j_set)
        for collection in relevant_pairings:
            path_collection_iterator = None
            # build up path_collection_iterator by taking the product of iterators
            # where each factor iterator is the ones that connect i to j
            ij_paths = (nx.all_simple_paths(
                self._underlying_graph, f"{cur_i}_0", f"{cur_j}_{self.sink_number[cur_j]}")
                for cur_i, cur_j in collection)
            path_collection_iterator = itertools.product(*ij_paths)
            # now path_collection is a tuple of k paths that fully connect i_set to j_set
            for path_collection in path_collection_iterator:
                seen_vertices: Set[Nat] = set({})
                is_vertex_disjoint = True
                my_paths = []
                for path in path_collection:
                    my_paths.append(list(path))
                    seen_vertices.isdisjoint(set(path))
                    if seen_vertices.isdisjoint(set(path)):
                        seen_vertices = seen_vertices.union(set(path))
                    else:
                        is_vertex_disjoint = False
                        break
                if is_vertex_disjoint:
                    all_path_weights = [self.path_weight(
                        path_nodes) for path_nodes in my_paths]
                    yield reduce(lambda x, y: x*y, all_path_weights, self.multiplicative_identity)

    def lindstrom_minor(self, i_set: Set[Nat], j_set: Set[Nat]) -> WeightType:
        """
        the minor where specify to include the rows and columns specified by i_set and j_set
        """
        minor_value = self.additive_identity
        had_contribution = False
        for contrib in self.vertex_disjoint_collection(i_set, j_set):
            had_contribution = True
            minor_value += contrib
        if not had_contribution:
            self.totally_connected = False
        return minor_value

    def draw(self,show_as_well=True) -> None:
        """
        draw with the edge weights
        """
        all_node_names = list(self._underlying_graph.nodes())
        all_positions = {z: self._underlying_graph.nodes[z]["position"]
                         for z in all_node_names}
        all_edge_labels = {z: str(self._underlying_graph.edges[z]["weight"])
                           for z in list(self._underlying_graph.edges())
                           if self._underlying_graph.edges[z]["weight"] != self.multiplicative_identity}
        nx.draw(
            self._underlying_graph,pos=all_positions,with_labels=False)
        nx.draw_networkx_edge_labels(
            self._underlying_graph,
            pos=all_positions,
            edge_labels=all_edge_labels)
        if show_as_well:
            plt.axis('off')
            plt.show()

    def totally_nonnegative(self, nonnegative_wt: Callable[[NT], bool]) -> bool:
        """
        is the associated weight matrix manifestly totally nonnegative
        by the constituent weights being all nonnegative
        """
        all_weights : Iterator[WeightType] = (self._underlying_graph[src][tgt]['weight']
                       for src, tgt in self._underlying_graph.edges())
        def real_nonnegative_wt(arg : Union[float,NT]) -> bool:
            """
            essentially same as nonnegative_wt
            """
            if isinstance(arg,(float,int)):
                return arg>=0
            return nonnegative_wt(arg)
        return all((real_nonnegative_wt(wt) for wt in all_weights))

    def positive(self, positive_wt: Callable[[NT], bool]) -> bool:
        """
        is the associated weight matrix manifestly totally positive
        by the constituent weights being all positive
        """
        assert self.totally_connected
        all_weights : Iterator[WeightType] = (self._underlying_graph[src][tgt]['weight']
                       for src, tgt in self._underlying_graph.edges())
        def real_positive_wt(arg : Union[float,NT]) -> bool:
            """
            essentially same as positive_wt
            """
            if isinstance(arg,(float,int)):
                return arg>=0
            return positive_wt(arg)
        return all((real_positive_wt(wt) for wt in all_weights))

U = TypeVar("U")

def pairings(i_set: Set[T], j_set: Set[U]) -> Iterator[Set[Tuple[T, U]]]:
    """
    generator that goes through all pairings of i_set with j_set
    """
    if len(i_set) != len(j_set):
        raise ValueError("Must have same number of elements")
    if len(i_set) == 0:
        yield set({})
        return
    xtemp = i_set.copy()
    i_popped = xtemp.pop()
    for my_i, my_j in {(i_popped, j) for j in j_set}:
        ytemp = j_set.copy()
        ytemp.remove(my_j)
        for remaining in pairings(xtemp, ytemp):
            remaining.add((my_i, my_j))
            yield remaining


@runtime_checkable
class Ordered(Protocol):
    """
    protocol for ordered
    """

    def __leq__(self, _other): ...


V = TypeVar('V', bound=Ordered)
W = TypeVar('W', bound=Ordered)


def pairings_nc(i_set: Set[V], j_set: Set[W]) -> Set[Tuple[V, W]]:
    """
    the unique pairing of i_set with j_set
    such that the minimum of i_set is matched with the minimum of j_set
    and so on all the way up
    """
    if len(i_set) != len(j_set):
        raise ValueError("Must have same number of elements")
    i_list = list(i_set)
    i_list.sort()
    j_list = list(j_set)
    j_list.sort()
    return set(zip(i_list, j_list))

# pylint:disable = too-many-arguments
def determinant(matrix : List[List[NT]], mul : NT,
                simplifier : Callable[[NT],NT],
                additive_identity : NT,
                negative_one : NT,
                simplify_freq : int=3,
                original_width : Optional[int] =None) -> NT:
    """
    determinant using only + and *
    so that it can work on complex object types
    that override __add__ and __mul__
    imported det functions assume primitive numeric type like float
    """
    width = len(matrix)
    if original_width is None:
        original_width = width
    if width == 1:
        to_return = mul * matrix[0][0]
        if original_width == width:
            to_return = simplifier(to_return)
        return to_return
    sign = negative_one
    answer = additive_identity
    for skip_col in range(width):
        m = []
        for jdx in range(1, width):
            buff = []
            for kdx in range(width):
                if kdx != skip_col:
                    buff.append(matrix[jdx][kdx])
            m.append(buff)
        sign *= negative_one
        answer = answer + mul * \
            determinant(m, sign * matrix[0][skip_col], simplifier,
                        additive_identity, negative_one, simplify_freq, original_width)
    if width % simplify_freq == original_width % simplify_freq:
        answer = simplifier(answer)
    return answer
