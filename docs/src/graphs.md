# Graphs

Currently only the square lattice is implemented

## Square lattice

Create a square lattice with periodic or open boundary conditions:

```@repl 1
	using LatticeModels
	lato = SQLattice_open{Int}(5,5)
	latp = SQLattice_periodic{Int}(5,5)
```

Index with integers, `CartesianIndex` or with a `Site` object.  The latter knows about boundary conditions and lattice topology.

```@repl 1
	lato .= reshape(1:length(lato),size(lato))
	latp .= reshape(1:length(lato),size(latp))
	lato[1,1] + latp[1,1]
	
	I = SQLattice_site(lato,1,1)
	J = SQLattice_site(latp,1,1)
	foreach_neighbour(I) do s println(s) end
	foreach_neighbour(J) do s println(s) end
```
 
