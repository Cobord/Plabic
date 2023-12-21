"""
A cluster consisting of a quiver with n
vertices
and an expression for each vertex
"""

from __future__ import annotations
from typing import Callable, Optional, TypeVar,Protocol,\
    Generic,List,Dict,Tuple,Set
from functools import reduce
import networkx as nx
import matplotlib.pyplot as plt
from sympy import Expr, Integer, symbols

T = TypeVar("T")

class Arithmetic(Protocol):
    """
    can be multiplied, added, divided and powered
    that is enough to do all the cluster
    coordinate changes
    """
    def __add__(self : T,other : T) -> T: ...
    def __mul__(self : T,other : T) -> T: ...
    def __pow__(self : T,power : int) -> T: ...
    def __truediv__(self : T,other : T) -> T: ...

NT = TypeVar("NT",bound=Arithmetic)
VertexLabel = str

def cluster_same(cluster1 : Dict[VertexLabel,NT],
                 cluster2: Dict[VertexLabel,NT],
                 force_simplify : Optional[Callable[[NT],NT]] = None) -> bool:
    """
    same cluster variables and same associated vertices in quiver
    """
    if cluster1.keys() != cluster2.keys():
        return False
    for which_vertex,v_unsubbed in cluster1.items():
        if not v_unsubbed != cluster2[which_vertex]:
            if force_simplify is not None:
                v_subbed = force_simplify(v_unsubbed)
                w_subbed = force_simplify(cluster2[which_vertex])
                if v_subbed == w_subbed:
                    continue
            return False
    return True

class Cluster(Generic[NT]):
    """
    a quiver and an instance of NT for each vertex
    """
    def __init__(self,
                 my_quiver : nx.MultiDiGraph,
                 initial_variables : List[NT],
                 multiplicative_identity : NT,
                 simplifier : Callable[[NT],NT]):
        if len(initial_variables) == 0:
            raise ValueError("There should be at least one generator")
        if my_quiver.number_of_nodes() != len(initial_variables):
            raise ValueError("There should be matching number of vertices and variables")
        first_vertex = list(my_quiver.nodes())[0]
        if not isinstance(first_vertex,VertexLabel):
            raise ValueError(
                f"The vertices should have type {VertexLabel}")
        for vertex in my_quiver.nodes():
            self_loop_count = my_quiver.number_of_edges(vertex,vertex)
            if self_loop_count>0:
                raise ValueError("No self loops")
        for (x_node,y_node) in my_quiver.edges():
            if my_quiver.number_of_edges(y_node,x_node)>0:
                raise ValueError("No oriented 2-cycles")
        self.my_quiver = my_quiver
        self.cluster : Dict[VertexLabel,NT] = {z : initial_variables[idx]
                        for idx,z in enumerate(my_quiver.nodes())}
        self.multiplicative_identity = multiplicative_identity
        self.simplifier = simplifier

    #pylint:disable=too-many-locals,too-many-branches
    def mutate(self,key : VertexLabel) -> None:
        """
        mutation at the vertex labelled key
        """
        if key not in self.cluster:
            raise ValueError(f"{key} was not a label for a vertex so we can't mutate there")
        mutating_variable = self.cluster[key]
        to_edges : Dict[VertexLabel,int] = {}
        from_edges : Dict[VertexLabel,int] = {}
        for (src_vertex,tgt_vertex) in self.my_quiver.edges(keys=False):
            if src_vertex==tgt_vertex:
                raise ValueError("There should be no self loops")
            if src_vertex==key:
                to_edges[tgt_vertex] = to_edges.get(tgt_vertex,0)+1
            if tgt_vertex==key:
                from_edges[src_vertex] = from_edges.get(src_vertex,0)+1
        summand1 = reduce(lambda acc,x:acc*(self.cluster[x[0]]**x[1]),
                          from_edges.items(),self.multiplicative_identity)
        summand2 = reduce(lambda acc,x:acc*(self.cluster[x[0]]**x[1]),
                          to_edges.items(),self.multiplicative_identity)
        self.cluster[key] = self.simplifier((summand1+summand2)/mutating_variable)

        xyz_modifications : Set[Tuple[VertexLabel,VertexLabel,VertexLabel]] = set()
        for (x_node,y_node) in self.my_quiver.in_edges(key):
            continuing_edges = ((y_node,z_node) for (y2_node,z_node)
                                in self.my_quiver.out_edges(key) if y2_node==y_node)
            for (_,z_node) in continuing_edges:
                xyz_modifications.add((x_node,y_node,z_node))
        for (x_node,y_node,z_node) in xyz_modifications:
            self.my_quiver.add_edge(x_node,z_node)
        to_flip : Set[Tuple[VertexLabel,VertexLabel,int]] = set()
        next_edge_index = 0
        for (src_vertex,tgt_vertex,edge_idx) in self.my_quiver.in_edges(key,keys=True):
            to_flip.add((src_vertex,tgt_vertex,edge_idx))
            next_edge_index = max(next_edge_index,edge_idx)
        for (src_vertex,tgt_vertex,edge_idx) in self.my_quiver.out_edges(key,keys=True):
            to_flip.add((src_vertex,tgt_vertex,edge_idx))
            next_edge_index = max(next_edge_index,edge_idx)
        for (src_vertex,tgt_vertex,edge_idx) in to_flip:
            self.my_quiver.remove_edge(src_vertex,tgt_vertex,edge_idx)
            next_edge_index += 1
            self.my_quiver.add_edge(tgt_vertex,src_vertex,next_edge_index)
        has_oriented_2_cycles : Set[Tuple[VertexLabel,VertexLabel]] = set()
        for (x_node,y_node) in self.my_quiver.edges():
            if self.my_quiver.number_of_edges(y_node,x_node)>0:
                if x_node>y_node:
                    x_node,y_node = y_node,x_node
                has_oriented_2_cycles.add((x_node,y_node))
        for (src_vertex,tgt_vertex) in has_oriented_2_cycles:
            src_2_tgt = self.my_quiver.number_of_edges(src_vertex,tgt_vertex)
            tgt_2_src = self.my_quiver.number_of_edges(tgt_vertex,src_vertex)
            net_src_2_tgt = src_2_tgt - tgt_2_src
            for _ in range(src_2_tgt):
                self.my_quiver.remove_edge(src_vertex, tgt_vertex)
            for _ in range(tgt_2_src):
                self.my_quiver.remove_edge(tgt_vertex, src_vertex)
            if net_src_2_tgt>=0:
                for _ in range(net_src_2_tgt):
                    self.my_quiver.add_edge(src_vertex, tgt_vertex)
            else:
                for _ in range(-net_src_2_tgt):
                    self.my_quiver.add_edge(tgt_vertex, src_vertex)

    def draw(self) -> None:
        """
        draw the quiver
        """
        nx.draw(self.my_quiver)
        plt.show()

    def my_antisymmetric(self) -> List[List[int]]:
        """
        the quiver as antisymmetric integer matrix
        """
        my_nodes = list(self.my_quiver.nodes())
        num_vertices = len(my_nodes)
        ret_val = [[0 for _ in range(num_vertices)] for _ in range(num_vertices)]
        for idx,src_node_name in enumerate(my_nodes):
            for jdx,tgt_node_name in enumerate(my_nodes):
                src_2_tgt = self.my_quiver.number_of_edges(src_node_name,tgt_node_name)
                tgt_2_src = self.my_quiver.number_of_edges(tgt_node_name,src_node_name)
                ret_val[idx][jdx] = src_2_tgt - tgt_2_src
        return ret_val

    def __eq__(self,other) -> bool:
        if not isinstance(other,Cluster):
            return False
        if list(self.my_quiver.nodes()) != list(other.my_quiver.nodes()):
            return False
        if not cluster_same(self.cluster,other.cluster):
            return False
        return self.my_antisymmetric() == other.my_antisymmetric()

    @staticmethod
    def make_type_an_cluster(my_n : int) -> Tuple[List[Expr],Cluster[Expr]]:
        """
        make the cluster for type An with sympy variables
        """
        if my_n<=0:
            raise TypeError("n must be a natural number greater than 0")
        variable_pool = ["a","b","c","d","x","y","z","w","v",
                         "f","g","h","l","m","n","o",
                         "p","q","r","s","t","u"]
        if my_n>len(variable_pool):
            raise ValueError(f"Only have {len(variable_pool)} variables available")
        variables_used = variable_pool[0:my_n]
        variable_list = list(symbols(",".join(variables_used),commutative=True))
        my_an_quiver = nx.MultiDiGraph()
        for cur_node in range(1,my_n+1):
            my_an_quiver.add_node(str(cur_node))
        for cur_node in range(1,my_n):
            my_an_quiver.add_edge(str(cur_node),str(cur_node+1))
        return_cluster = Cluster[Expr](my_an_quiver,variable_list,Integer(1),lambda z: z.simplify())
        return variable_list, return_cluster
