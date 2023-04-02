# Sandpiles.jl -- Sandpile models
#
# Copyright (C) 2023 Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For details, see the file LICENSE in the root directory, or
# check <https://www.gnu.org/licenses/>.

using Random

export BoundaryType, b_open, b_semiopen
@enum BoundaryType b_open b_semiopen

abstract type Boundary_conditions end
struct Open_boundaries <: Boundary_conditions end
struct Semiopen_boundaries <: Boundary_conditions end

export BTW_sandpile, run!

mutable struct BTW_sandpile{RNG_T<:AbstractRNG, BC_T<:Boundary_conditions}
    L::Int
    z::Matrix{Int}   # A boundary layer is used to implement boundary conditions
    zc::Int
    ztot::BigInt
    time::Int
    BC::BC_T
    RNG::RNG_T
    z_N             # These will hold views to the boundary layers
    z_S
    z_E
    z_W
    z_in            # View to the actual sandpile heights
end

export Sandpile_history
"""
    struct Sandpile_history

struct to hold historical information of the sandpile: average height
and avalanche characteristics (start, duration, size and energy).
Will be properly updated if passed to the `run` method.
"""
mutable struct Sandpile_history
    time::Vector{Int}
    zav::Vector{Float64}
    avalanche_start::Vector{Int}
    avalanche_duration::Vector{Int}
    avalanche_size::Vector{Float64}
    avalanche_energy::Vector{Int}
end

"Return an empty history object"
Sandpile_history()=Sandpile_history(Int[],Float64[],Int[],Int[],Float64[],Int[])

"""
Return a history object initialized with the average
height of the given sandpile, and empty avalanche
information
"""
Sandpile_history(pile::BTW_sandpile) =
    Sandpile_history([pile.time],[pile.ztot/pile.L^2],Int[],Int[],Float64[],Int[])

"""
    BTW_sandpile(L,zc::Int;seed=443,boundaries::BoundaryType=open)

Create a 2-d BTW a new sandpile with random heights, setting size `L`
and the critical height `zc`.  'boundaries' can be 'open' or
'semiopen'.  Optionally, a seed for the number generator can be
passed.

To simulate the sandpile, create a BTW_sandpile object and optionally
a Sandpile_history, and call `run!`.

The implementation of the BTW sandpile model, roughly follows

H. J. Jensen, _Self-Organized Criticality_, Cambridge University Press,
Cambridge (1998)
"""
function BTW_sandpile(L,zc::Int;seed=443,boundaries::BoundaryType=open)
    RNG = Random.Xoshiro(seed)
    z=zeros(Int,L+2,L+2)
    z_N=@view z[1,:]
    z_S=@view z[L+2,:]
    z_E=@view z[:,1]
    z_W=@view z[:,L+2]
    z_in=@view z[2:L+1,2:L+1]
    s = boundaries==b_open ?
        BTW_sandpile(L,z,zc,BigInt(0),0,Open_boundaries(),RNG,z_N,z_S,z_E,z_W,z_in) :
        BTW_sandpile(L,z,zc,BigInt(0),0,Semiopen_boundaries(),RNG,z_N,z_S,z_E,z_W,z_in)
    rand!(RNG,s.z_in,1:zc)
    s.ztot=sum(s.z_in)
    return s
end

"""
    run!(pile::BTW_sandpile;
                  steps::Int=1,history::Union{Sandpile_history,Nothing}=nothing)

Do at least `steps` steps of dynamical evolution for the sandpile
`pile`.  If at the requested number of steps an avalanche
is ongoing, continue until avalanche is finished and all
sites are subcritical.
"""
function run!(pile::BTW_sandpile;
              steps::Int=1,history::Union{Sandpile_history,Nothing}=nothing,
              reference_conf=nothing)

    step0=pile.time
    while pile.time<step0+steps
        ir=rand(pile.RNG,CartesianIndices(pile.z_in))
        pile.z_in[ir]+=1
        pile.ztot+=1
        if pile.z_in[ir]>pile.zc avalanche(pile,ir,history)
        else
            pile.time+=1
            if !isnothing(history)
                push!(history.time,pile.time)
                push!(history.zav,pile.ztot/pile.L^2)
                if !isnothing(reference_conf)
                    compute_overlap(pile,reference_conf,history)
                end
            end
        end
    end

    @assert critical_sites(pile) == 0
    
end

critical_sites(pile::BTW_sandpile) = sum(pile.z .> pile.zc)

Siteset = Set{Tuple{Int,Int}}

"""
Develop an avalanche until no supercritial sites remain.  This
is called from `run!` when some site becomes supercritical.
"""
function avalanche(pile::BTW_sandpile,scsite::CartesianIndex{2},
                   history::Union{Sandpile_history,Nothing})

    if !isnothing(history) push!(history.avalanche_start,pile.time) end
    supercritical_sites::Siteset=Set([ parentindices(@view(pile.z_in[scsite])) ])
    avener=0
    cluster=Set(supercritical_sites)
    while length(supercritical_sites)>0
        avener+=length(supercritical_sites)
        relax(pile,supercritical_sites,pile.BC)
        pile.time+=1
        if !isnothing(history)
            push!(history.time,pile.time)
            push!(history.zav,pile.ztot/pile.L^2)
        end
        update_sc!(pile,supercritical_sites)
        union!(cluster,supercritical_sites)
    end
    if !isnothing(history)
        push!(history.avalanche_duration,pile.time-history.avalanche_start[end])
        push!(history.avalanche_energy,avener)
        push!(history.avalanche_size,csize(cluster))
    end

end

"Compute size of cluster"
function csize(cluster::Set)
    CM=Float64[0, 0]
    N=length(cluster)
    for I in cluster CM+=[ I[1], I[2] ] end
    CM/=N
    l=0.
    for I in cluster l+= sqrt( (I[1]-CM[1])^2 + (I[2]-CM[2])^2 ) end
    l/=N
    return l
end

"""
Relax all the sites in the `supercritical_sites` collection using semi-open
boundary conditions.  Sites are not checked again, it is assumend that
`supercritical_sites` holds the right list.  May generate
new supercritical sites.
"""
function relax(pile::BTW_sandpile,supercritical_sites,_::Semiopen_boundaries)
    for (jx,jy)  in supercritical_sites
            pile.z[jx  , jy] -= 4
            pile.z[jx+1, jy] += 1
            pile.z[jx-1, jy] += 1
            pile.z[jx, jy+1] += 1
            pile.z[jx, jy-1] += 1
    end

    # Bounce back grains from the closed sides
    pile.z[2,:] .+= pile.z_N
    pile.z[:,2] .+= pile.z_E
    pile.z_N .= 0
    pile.z_E .= 0

    # Drop grains from the other boundaries
    pile.ztot-=sum(pile.z_S)
    pile.ztot-=sum(pile.z_W)
    pile.z_S .= 0
    pile.z_W .= 0
end

"""
Relax all the sites in the `supercritical_sites` collection using open
boundary conditions.  Sites are not checked again, it is assumend that
`supercritical_sites` holds the right list.  May generate
new supercritical sites.
"""
function relax(pile::BTW_sandpile,supercritical_sites,_::Open_boundaries)
    for (jx,jy)  in supercritical_sites
            pile.z[jx  , jy] -= 4
            pile.z[jx+1, jy] += 1
            pile.z[jx-1, jy] += 1
            pile.z[jx, jy+1] += 1
            pile.z[jx, jy-1] += 1
    end

    # Drop grains over boundaries
    pile.ztot -= sum(pile.z_N)
    pile.ztot -= sum(pile.z_S)
    pile.ztot -= sum(pile.z_E)
    pile.ztot -= sum(pile.z_W)

    pile.z_N .= 0
    pile.z_S .= 0
    pile.z_E .= 0
    pile.z_W .= 0

end

@inline check(pile,sc,ix,iy) =  if pile.z[ix,iy]>pile.zc union!(sc,[(ix,iy)]) end

function update_sc!(pile::BTW_sandpile,supercritical_sites)
    sc=copy(supercritical_sites)
    empty!(supercritical_sites)
    for (ix,iy) in sc
        check(pile,supercritical_sites,ix,iy)
        if ix>1 check(pile,supercritical_sites,ix-1,iy) end
        if ix<=pile.L check(pile,supercritical_sites,ix+1,iy) end
        if iy>1 check(pile,supercritical_sites,ix,iy-1) end
        if iy<=pile.L check(pile,supercritical_sites,ix,iy+1) end
    end
end
