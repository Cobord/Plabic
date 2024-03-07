# Plabic
Cluster algebraic/geometric structures related to plabic graphs

## Requirements

matplotlib, networkx, numpy, sympy


# Framed 2 Disks

```
from plabic import FramedDiskConfig
e = FramedDiskConfig([((0,.5),.2,0),((-.5,0),.3,3.14)],outer_circle=None)
```

This means there is an internal circle at (0,.5) with radius .2 and the framing is specified by the angle 0 relative to the positive x-axis.

There is also another internal circle at (-.5,0) with radius .3 and the framing is specified by the angle 3.14 relative to the positive x-axis.

The external circle is optional and defaults to the standard unit circle around (0,0) with radius 1. You can specify it with the same
format as center,radius,angle.

This can be drawn and undergo operadic substituion, but those are internal implementaion details for drawing and doing operadic substitutions on PlabicGraph.

# Plabic Diagram

## Plabic Graph Builder

```
from plabic import PlabicGraphBuilder
builder = PlabicGraphBuilder()
builder.set_num_external(5)
builder.set_internal_circles_nums([1,2])
builder.add_external_bdry_vertex("ext1",0,"int1",extras=None)
builder.add_internal_bdry_vertex("bdry1",1,1,"int1",extras=None)
builder.add_internal_vertex("int1",BiColor.RED,["ext1","int2","bdry1"],extras=None)
....
builder.set_circles_config(e)
p = builder.build()
```

Here we build a Plabic Graph by first saying there are 5 vertices on an external boundary.

There are then stated to be 2 internal boundaries with 1 and 2 vertices on them respectively.

Then we add that the 0th vertex on the external boundary is named "ext1" and connects to "int1".

After that we add that the 1st vertex on the 1st internal boundary is named "bdry1" and connects to "int1".

We then add the internal vertex named "int1" which is colored RED and connects to "ext1","int2","bdry1" in clockwise order.

Then there are more steps including telling the information about "int2"'s color and connections.

We then state that the boundaries are the circles given as above. This does not really hold in this case bc we did not
give positions for all the vertices, but we could have done that. They would be in the extras which we put as None above.
Those extras are passed as a dictionary like 

```
extras = {"position": (.5,.3),"my_perfect_edge":0,...}
```

After all that we can build and all this information goes into the construction of a PlabicGraph. We could also pass the
information to PlabicGraph construction immediately, but that is beyond the scope of most common use.

There is also a method for set_multi_edges but that is only relevant when there are multiple edges between two vertices and/or self loops.

You are also allowed to rename vertices as long as the new name would not cause a name collision with a name already present.

## Plabic Graph

```
from plabic import PlabicGraph
p.remove_prop("my_perfect_edge")
_ = p.is_trivalent("int2")
_ = p.one_color_trivalent(BiColor.RED)
_ = p.boundary_connects_to_color(BiColor.GREEN)
_ = p.is_bipartite()
_ = p.nodes_connected("int1","int2")
_ = p.num_connecting_edges("int1","int2")
_ = p.get_color("int1")
_ = p.opposite_colors("int1","int2")
lollipop_color,destination = p.bdry_to_bdry("ext1")
p.draw()
```

1. Removing the property my_perfect_edge from all the vertices
2. Querying if the vertex int2 is trivalent
3. Querying if all RED vertices are trivalent
4. Querying if boundary vertices only connect to GREEN internal vertices
5. Querying if the graph is bipartite of RED and GREEN
6. Querying if there is any edge connecting int1 and int2
7. Querying how many edges connect int1 and int2
8. Querying the color of the vertex named int1 and boundary vertices have no color
9. Querying whether int1 and int2 are both colored and that they are not the same color
10. Perform the decorated permutation starting at ext1 and give the destination and the decoration. There is a helper for doing one step and it is called follow_rules_of_road
11. Draw the graph with matplotlib,
      There are many options of what colors to use and how to draw the edges and boundaries, but one can use the defaults.

### Moves on Plabic Graphs

- greedy_shrink does M2 and M3 moves as much as possible, if the vertices have positions
    then the M3 moves don't occur because we don't know where to place the collapsed vertex
    they have to be done separately with more information provided
- square_move is given 4 vertex names. If those form a square of trivalent internal vertices
    and they have alternating colors, swap the colors around. If any of these preconditions fail,
    then the output is False and an explanation along with no mutation of the Plabic Graph
- flip_move is given 2 vertex names. If those form an edge of trivalent internal vertices
    and they have the same colors, reconnect the 4 remaining half edges in the other way.
    We can also instruct how to transform the extra data on the two vertices that are changing.
    The extra data on the two new vertices like precise position are given by the output of a provided function.
    This function is optional and defaults to not having any extra data on the new vertices. In particular, there is no longer a specific position on the plane for each vertex in that case.
    If any of these preconditions fail, then the output is False and an explanation along
    with no mutation of the Plabic Graph
- remove_bivalent takes in a vertex name. If that is a bivalent vertex, then remove it. If any of these preconditions fail,
    then the output is False and an explanation along with no mutation of the Plabic Graph
- insert_bivalent receives 2 vertex names, a color, a name and a function that produces the extra data of the new vertex.
    If the two vertices provided have a single edge between them, then a new vertex with the specified name and color is inserted onto that edge.
    The extra data like it's precise position is given by the output of the provided function. The function is optional and defaults to not having any
    extra data on that new vertex. In particular, there is no longer a specific position on the plane for each vertex.
    If any of these preconditions fail, then the output is False and an explanation along
    with no mutation of the Plabic Graph
- contract_edge receives 2 vertex names, a name and a function that produces the extra data of the new vertex.
    If the two vertices provided have a single edge between them and they are the same color, then that edge collapses into a single new vertex.
    That new vertex has the specified name and color. The extra data like it's precise position is given by the output of the provided function.
    The function is optional and defaults to not having any extra data on that new vertex. In particular, there is no longer a specific position on the plane for each vertex.
    If any of these preconditions fail, then the output is False and an explanation along
    with no mutation of the Plabic Graph        
- split_vertex receives a vertex name, a pair of integers, two new vertex names and a function that produces the extra data of the new vertices.
    the half edges out of the given vertex are split into two groups. For example, if the pair was (1,5) then 1 through 5 would go to one vertex and the others would go to the other.
    These indices are treated periodically, so we could also do (5,1) and get a different splitting. These become new vertices with the provided new names.
    The extra data for the new vertices is determined by that function and again that function is optional. Again if any of the preconditions fail,
    we get a False output and an explanation along with no mutation to the Plabic Graph.
- operad_compose takes two Plabic Graphs, and an index which_internal_disk and glues the second in
    at the appropriate internal boundary of the former. If any of the conditions for a well defined gluing fail,
    then the output is False and an explanation along with no mutation of the Plabic Graph. If there are circle configurations as in FramedDiskConfig, then
    methods described there are used and the updated FramedDiskConfig for the glued configuration is set.

# Le Diagram

```
from plabic import LeDiagram
my_Le = LeDiagram([[0,1,0,1,0],[1,1,0,1],[0,0],[0,1]])
second_col = my_Le.column_height(1)
for (a,b) in my_Le.nw_path((3,1),False):
    pass
my_grassman_necklace = my_Le.to_grassmann_necklace(7,7)
my_Le_again = LeDiagram.from_grassmann_necklace(my_grassmann_necklace)
p = my_Le.to_plabic()
_ = my_Le.to_special_bruhat_interval(bounding_k,bounding_n)
```

One may also construct the Le Diagram with True/False filling instead of 1/0.

We can query the heights of the columns. The lengths of the rows are just len(my_Le.filling[row_number])

We can iterate through the path that always going strictly northwest which starts at the nearest 1 strictly/weakly northwest of a starting cur_loc.
That is helped by a nw_path_next method.

Conversion to Grassmann necklaces proceeds by providing the bounding rectangle (first k second n). The output Grassmann necklace is a list of sets of integers. Each neighboring pair of sets satisfies the condition for being a Grassmann necklace which indicates precisely how those two sets can differ.

Conversely we can go back from a Grassmann necklace to a Le Diagram. One can optionally provide the n and k which are the length of the necklace and the sizes of the sets (they are all the same size by the Grassmann condition).

The Plabic Graph produced from the Le Diagram follows the standard rules based on positions of 1's. The vertices are positioned inside the boxes of the Young diagram as they would appear in English notation with all the boxes being unit squares.

There is also a method which produces a pair of permutations in interval format. If we call them v and w respectively, then v<=w and so we put them
in a Bruhat Interval. The w is k-Grassmannian for the k which is specified by the bounding box. They are both permutations of n which is also specified by the bounding box as the sum of the width and height.

# Triangulation

```
from plabic import Triangulation
from math import sin, cos, pi as PI
N = 8
my_diagonals = [(4,7),(2,8),(2,4),(2,7),(4,6)]
RADIUS = 2
t = Triangulation([ (RADIUS*cos(2*PI*which/N),RADIUS*sin(2*PI*which/N))
                       for which in range(N)], [(x-1,y-1) for x,y in my_diagonals])
change, new_diag = t.quad_flip((1,3))
p = t.to_plabic()
```

We create a triangulation of a convex N-gon by specfying the locations of the vertices 1 through N and the diagonals. Here the diagonal (4,7) means vertex 4 and 7 are connected by a diagonal (with the -1's for zero indexing). From the diagonals the triangles are all deterrmined as well as all the respective quadrilaterals each diagonal is part off with it's two neighboring triangles.

Then we can do diagonal flip moves on the triangulation by specifying which diagonal to flip. Here we are flipping the diagonal which connects vertex 2 and vertex 4 (off by 1 due to zero indexing). If the specified input was not a diagonal of a quadrilateral that could be flipped, then the first output indicating whether the Triangulation has changed will be False. Otherwise we will get True and the new diagonal which was the other diagonal of the relevant quadrilateral.

This produces a special Plabic Graph. There are internal and external vertices for each vertex of the polygon and internal vertices for each triangle. The external vertices corresponding to the polygon vertices are connected to the corresponding internal vertex. The vertices for the triangles are connected to the vertices for their corners. This Plabic Graph gets a perfect orientation. If all the vertices of the N-gon were on a circle like above, then a FramedDiskConfig with 0 internal circles is also produced. When the Plabic Graph produced from the Triangulation is drawn, that circle will also be drawn passing through all the boundary vertices.

# Double Wiring Diagrams

```
from plabic import WiringDiagram
w = WiringDiagram(5, [1,-2,3,-3])
for (c1,c2) in w.chamber_minors():
  pass
p = w.to_plabic()
```

for creating a Double wiring diagram for S_5 x S_5 which consists of s_1 s_{-2} s_3 s_{-3} where negative indices mean using the other strands of the DOUBLE wiring diagram

c1,c2 in the chamber minors are sets of integers. They indicate taking rows that are in c1 and columns that are in c2. For how they are determined by a double wiring diagram, see the original Fomin-Zelivinsky papers.

to_plabic gives the corresponding Plabic graph. The vertices are automatically given positions in the plane and if the diagram ends up being bipartite, then it the Plabic graph produced also includes a pre-chosen perfect orientation.

# Planar Diagram

```
from plabic import PlanarNetwork
p = PlanarNetwork(3, edge_list=[[(3, 2, A), (3, 3, ONE)], [(3, 2, C), (2, 1, B)], [
                      (2, 2, E), (1, 1, D)], [(3, 3, F), (2, 3, G), (1, 2, H)], [(2, 3, I)]],
                      multiplicative_identity=ONE,
                      additive_identity=ZERO)
a_12 = p.weight_matrix(1,2)
Delta_12_23 = p.lindstrom_minor(set([1,2]),set([2,3]))
assert p.totally_nonnegative(lambda letter : is_nonnegative_variable(letter))
assert p.positive(lambda letter : is_positive_variable(letter))
```

Consider the example in https://arxiv.org/pdf/math/9912128.pdf figure 1. We can see there are 3 horizontal lines as the first input states. The second input is the edge_list which describes the edges as read left to right in the figure along with variables A-I for their weights, and some being labelled with weight ONE. We can specify the multiplicative and additive identity. They default to 1.0 and 0.0 but if we want to work symbolically or with other number systems this provides the capability to do so. The weights simply need to implement the protocol of being able to be added and multiplied.

Alternatively instead of specifying edge_list one can also give a chip_word and chip_weights. The former is a list of tuples with the first part being ChipType.UP,DOWN,FLAT and the second part being a natural number to indicates for which pair of neighboring wires this "chip" is relevant. The latter is a list of the same length with the weights on the respective chips. In this construction, the diagram is more spaced out even if multiple chips can be placed next to each other vertically.

We can then determine the ij entry of the weight matrix which is given by a sum of products expression in the edge weights for paths connecting i on one side of the diagram to j on the other.

Lindstrom Minor uses the Lindstrom lemma to compute the specified minor of the weight matrix via a similar sum of products formula with systems of paths connecting the vertices in one set on the left side to vertices in the other set on the right side.

The 2 totally_nonnegative/positive query methods take the diagram and a function which says whether or not a given edge weight is nonnegative/positive and deduces whether the weight matrix also is totally nonnegative/totally positive. This is when we have defined A-I via sympy to be unspecified variables that we will later interpret to only take nonnegative/positive values.

When calling the draw method a figure is produced with the diagram. The edges are labelled with their weights (except if the weight is the multiplicative identity which is treated as default and doesn't need to be explicitly shown).

# (Bounded) Affine Permutations

There are many ways to construct these.

- We can specify where several integers go enough to completely pin down everything with that and shift equivariance.
- We may use a list of the Coxeter generators where the Coxeter graph is affine An
- We may use juggling patterns
- We may specify a decorated trip permutation and produce a bounded affine permutation from that
- We may multiply, invert and concatenate to build from smaller examples

There are several generators which yield all sorts of (bounded) affine permutations with various properties

- All bounded affine permutations for a given n
- Everything in S_n up to given Coxeter lengths (inessentially embedded into the affine case for common interface)
- Everything in S_k times S_{n-k} (inessentially embedded into the affine case for common interface)
- Everything in the affine symmetric group up to a given Coxeter length
- Iterating through Q_{k,N} by lengths of the two factors (which form endpoints of a Bruhat interval)
- Shifts of such Bruhat intervals using the other generators

The to_plabic method uses the BCFW bridge construction to construct a plabic graph from a bounded affine permutation.

# Cluster

```
from plabic import Cluster,cluster_same
A,B,C = symbols("a,b,c",commutative=True)
my_a3_quiver = nx.MultiDiGraph()
my_a3_quiver.add_node("1")
my_a3_quiver.add_node("2")
my_a3_quiver.add_node("3")
my_a3_quiver.add_edge("1","2")
my_a3_quiver.add_edge("2","3")
c = Cluster[Expr](my_a3_quiver,[A,B,C],multiplicative_identity=Integer(1),simplifier=lambda z: z.simplify())
c.mutate("1")
assert cluster_same(c.cluster,{"1":(B+1)/A,"2":B,"3":C})
c.draw()
```

Create a cluster by specifying a quiver and expressions associated to each vertex. One also specifies a multiplicative identity which has the same generic type (Here all the variables and multiplicative identity are Expr from sympy).
We also pass a simplifier so that when mutations happen the repeated divisions are simplified so that we get back to a nice rational expression.

In order to mutate specify the name of the vertex you wish to mutate at. The object changes my_quiver and cluster variables (simply called cluster) appropriately. The cluster variables are given as a dictionary
so we can easily pick out the one corresponding to each vertex of the quiver.

cluster_same compares the two clusters and can also be given an additional function which can do more simplification. This will be needed only if the simplifier in the Cluster initialization was too weak. 

The draw method simply draws the underlying quiver.

# TNN Grassmannian

MinorSignConstraints allows one to specify which minors of a matrix are 0, positive, negative, etc.
The minors are specified as which columns are selected from the matrix. The total number of columns is not fixed,
just that it has to be at least some value for all the constraints to make sense.
You can change the number of columns in your example, and the additional columns will be unconstrained.

TNNGrassChart is where you actually work with the symbolic matrix representing a R_>0^d chart in Mat(k,n)
You specify the names of the variables which are all treated as ranging over positive reals
The constraints are imposed and we can enforce that all the sign constraints on minors are manifestly satisfied 

# BCFW Cells

[As described in this paper](https://arxiv.org/pdf/2402.15568.pdf "A cluster of results on amplitudehedron tiles")
Uses the operad composition for the BCFW product operation.