# Plabic
Cluster algebraic/geometric structures related to plabic graphs

# Requirements

matplotlib
networkx
numpy
sympy

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
