# Plabic
Cluster algebraic/geometric structures related to plabic graphs

# Requirements

matplotlib
networkx
numpy
sympy

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

# Triangulation

# Planar Diagram

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

# 
