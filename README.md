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

This can be drawn and undergo operadic substituion, but that is internal implementaion detail for drawing and doing operadic substitutions on PlabicGraph.

# Plabic Diagram

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
