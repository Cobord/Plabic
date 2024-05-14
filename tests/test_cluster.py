"""
test for clusters
"""
from typing import Dict
from sympy import Expr,symbols,Integer
import networkx as nx
import numpy as np
from plabic import Cluster,cluster_same

def my_cluster_same(cluster1 : Dict[str,Expr],
                 cluster2: Dict[str,Expr],
                 one_symbol=symbols("one",commutative=True),
                 zero_symbol=symbols("zero",commutative=True)) -> bool:
    """
    same cluster variables and same associated vertices in quiver
    """
    return cluster_same(cluster1,cluster2,lambda z: z.subs({one_symbol : 1.0,zero_symbol : 0.0}))

def test_a2():
    """
    a2 quiver
    """
    #pylint:disable=invalid-name
    A,B = symbols("a,b",commutative=True)
    my_a2_quiver = nx.MultiDiGraph()
    my_a2_quiver.add_node("1")
    my_a2_quiver.add_node("2")
    my_a2_quiver.add_edge("1","2")
    c = Cluster[Expr](my_a2_quiver,[A,B],Integer(1),lambda z: z.simplify())
    assert np.array_equal(c.my_antisymmetric(),[[0,1],[-1,0]])
    assert my_cluster_same(c.cluster,{"1":A,"2":B})
    c.mutate("1")
    assert np.array_equal(c.my_antisymmetric(),[[0,-1],[1,0]])
    assert my_cluster_same(c.cluster,{"1":(B+1)/A,"2":B})
    c.mutate("1")
    assert np.array_equal(c.my_antisymmetric(),[[0,1],[-1,0]])
    assert my_cluster_same(c.cluster,{"1":A,"2":B})
    c.mutate("1")
    c.mutate("2")
    assert np.array_equal(c.my_antisymmetric(),[[0,1],[-1,0]])
    assert my_cluster_same(c.cluster,{"1":(B+1)/A,"2":(1+A+B)/(A*B)})
    c.mutate("1")
    c.mutate("2")
    c.mutate("1")
    assert np.array_equal(c.my_antisymmetric(),[[0,-1],[1,0]])
    assert my_cluster_same(c.cluster,{"1":B,"2":A})

def test_make_a2():
    """
    a2 quiver with static method
    """
    #pylint:disable=invalid-name
    ([A,B],c) = Cluster.make_type_an_cluster(2)
    assert np.array_equal(c.my_antisymmetric(),[[0,1],[-1,0]])
    assert my_cluster_same(c.cluster,{"1":A,"2":B})
    c.mutate("1")
    assert np.array_equal(c.my_antisymmetric(),[[0,-1],[1,0]])
    assert my_cluster_same(c.cluster,{"1":(B+1)/A,"2":B})
    c.mutate("1")
    assert np.array_equal(c.my_antisymmetric(),[[0,1],[-1,0]])
    assert my_cluster_same(c.cluster,{"1":A,"2":B})
    c.mutate("1")
    c.mutate("2")
    assert np.array_equal(c.my_antisymmetric(),[[0,1],[-1,0]])
    assert my_cluster_same(c.cluster,{"1":(B+1)/A,"2":(1+A+B)/(A*B)})
    c.mutate("1")
    c.mutate("2")
    c.mutate("1")
    assert np.array_equal(c.my_antisymmetric(),[[0,-1],[1,0]])
    assert my_cluster_same(c.cluster,{"1":B,"2":A})

def test_square():
    """
    affine a3 quiver
    """
    #pylint:disable=invalid-name
    A,B,C,D = symbols("a,b,c,d",commutative=True)
    my_aa3_quiver = nx.MultiDiGraph()
    my_aa3_quiver.add_node("1")
    my_aa3_quiver.add_node("2")
    my_aa3_quiver.add_node("3")
    my_aa3_quiver.add_node("4")
    my_aa3_quiver.add_edge("1","2")
    my_aa3_quiver.add_edge("2","3")
    my_aa3_quiver.add_edge("3","4")
    my_aa3_quiver.add_edge("4","1")
    c = Cluster[Expr](my_aa3_quiver,[A,B,C,D],Integer(1),lambda z: z.simplify())
    c.mutate("1")
    assert c.my_quiver.number_of_edges() == 5
    assert c.my_quiver.number_of_edges("2","1")==1
    assert c.my_quiver.number_of_edges("1","4")==1
    assert c.my_quiver.number_of_edges("2","3")==1
    assert c.my_quiver.number_of_edges("3","4")==1
    assert c.my_quiver.number_of_edges("4","2")==1
    c.mutate("1")
    assert c.my_quiver.number_of_edges() == 4
    assert c.my_quiver.number_of_edges("1","2")==1
    assert c.my_quiver.number_of_edges("2","3")==1
    assert c.my_quiver.number_of_edges("3","4")==1
    assert c.my_quiver.number_of_edges("4","1")==1

def test_a3():
    """
    a3 quiver
    """
    #pylint:disable=invalid-name
    A,B,C = symbols("a,b,c",commutative=True)
    my_a3_quiver = nx.MultiDiGraph()
    my_a3_quiver.add_node("1")
    my_a3_quiver.add_node("2")
    my_a3_quiver.add_node("3")
    my_a3_quiver.add_edge("1","2")
    my_a3_quiver.add_edge("2","3")
    c = Cluster[Expr](my_a3_quiver,[A,B,C],Integer(1),lambda z: z.simplify())
    assert my_cluster_same(c.cluster,{"1":A,"2":B,"3":C})
    c.mutate("1")
    assert my_cluster_same(c.cluster,{"1":(B+1)/A,"2":B,"3":C})
    c.mutate("1")
    assert my_cluster_same(c.cluster,{"1":A,"2":B,"3":C})
    c.mutate("1")
    c.mutate("2")
    assert my_cluster_same(c.cluster,{"1":(B+1)/A,"2":(A+(1+B)*C)/(A*B),"3":C})
    c.mutate("3")
    assert my_cluster_same(c.cluster,{"1":(B+1)/A,"2":(A+(1+B)*C)/(A*B),"3":((1+B)*(A+C))/(A*B*C)})
