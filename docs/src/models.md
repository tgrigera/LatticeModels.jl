# Models

## Ising model

Ising model implemented over a generic graph.  For example, for the square lattice

```@example 1
using LatticeModels

IS = Ising(SQLattice_periodic,20,20;T=2.7)
println("Energy $(IS.E), magnetization $(IS.M)")
```

The single-spin Metropolis and Wolff cluster Monte Carlo algorithms are implemented:

```@example 1
E,M = Metropolis!(IS,steps=100,save_interval=10)
E,M = Wolff!(IS,steps=100,save_interval=10)
```

### API

```@docs
Ising
Metropolis!
Wolff!
```


### Analytical results

The following functions return exact values (in the thermodynamic limit), using Onsager's results.  These are taken from the [IsingModels.jl](https://github.com/cossio/IsingModels.jl) package by [
Jorge Fernandez-de-Cossio-Diaz](https://github.com/cossio).

```@docs
Ising_SQ_critical_temperature
Onsager_magnetization
Onsager_internal_energy
Onsager_heat_capacity
```

## Gaussian model
