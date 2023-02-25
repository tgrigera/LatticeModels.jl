# IsingModel.jl -- The Ising model on a generic graph
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
import SpecialFunctions

mutable struct Ising{Graph,RNG_T<:AbstractRNG}
    h::Float64
    T::Float64
    steps::Int64

    β::Float64
    βh::Float64
    M::Int64        # Magnetisation
    E0::Int64       # Energy w/o field
    E::Float64      # Actual energy

    conf::Graph
    probt::Matrix{Float64}
    RNG::RNG_T
end

function Ising(graph;seed=33,T=1.,h=0.)
    #RNG=Random.MersenneTwister(seed)
    RNG=Random.Xoshiro(seed)
    z=connectivity(graph)
    probt=zeros(Float64,z+1,2)
    I=Ising(h,0.,0,0.,0.,0,0,0.,graph,probt,RNG)
    set_temperature!(I,T)
    rand!(I.conf.nodes,[-1,1])
    set_energy_mag!(I)
    return I
end

function set_temperature!(IS::Ising,T)
    IS.T=T
    IS.β=1. /T
    IS.βh=IS.β*IS.h
    z=connectivity(IS.conf)
    for S=0:1
      for H=0:z
	DE=( 2*(-z+2*H) + 2*IS.h ) * (2*S-1);
	IS.probt[H+1,S+1] = DE<0  ?  2 : exp(-IS.β*DE)
      end
    end
end

"""
    βc

Critical temperature of the 2-dimensional infinite lattice Ising model determined
by Onsager.
"""
const βc = log1p(√2) / 2

"""
    onsager_magnetization(β)

Analytical magnetization of the two-dimensional Ising model
found by Onsager in the thermodynamic limit.
"""
function onsager_magnetization(β::Real)
    @assert β ≥ 0
    return max(1 - csch(2β)^4, 0)^(1/oftype(β, 8))
end

"""
    onsager_internal_energy(β)

Internal energy per site of the two-dimensional Ising model, derived by
Onsager in the thermodynamic limit.
"""
function onsager_internal_energy(β::Real)
    k = 2tanh(2β) / cosh(2β)
    j = 2tanh(2β)^2 - 1
    K = SpecialFunctions.ellipk(k^2)
    return -coth(2β) * (1 + 2/π * j * K)
end

"""
    onsager_heat_capacity(β)

Specific heat capacity of the two-dimensional Ising model, derived by
Onsager in the thermodynamic limit.
"""
function onsager_heat_capacity(β::Real)
    k = 2tanh(2β) / cosh(2β)
    K = SpecialFunctions.ellipk(k^2)
    E = SpecialFunctions.ellipe(k^2)
    j = 2tanh(2β)^2 - 1
    return β^2 * coth(2β)^2 * (2/π) * (((j - 1//2)^2 + 7//4) * K - 2E - (1 - j) * π / 2)
end

magnetization(IS::Ising{T}) where T = sum(IS.conf.nodes)

function energy_h0(IS::Ising{T}) where T
    E=Ref(0)
    foreach_bond(IS.conf) do a,b
        E[]+=a*b
    end
    return E[]
end

function set_energy_mag!(IS::Ising{T}) where T
    IS.E0=energy_h0(IS)
    IS.M=magnetization(IS)
    IS.E= IS.E0 - IS.h*IS.M
end

function Metropolis_sweep!(IS::Ising)
    site=random_site(IS.conf)
    localH=Ref{Int8}(0)
    for _=1:length(IS.conf)
        random_site!(IS.RNG,site)
        localH[]=0
        foreach_neighbour(site) do S
            localH[]+=S
        end
        S=IS.conf[site]
        p=IS.probt[(localH[] + connectivity(IS.conf))÷2+1,(S+1)÷2+1]
        if p>1. || rand(IS.RNG)<p

            IS.E0+=2*localH[]*S
            IS.M-=2*S
            IS.conf[site]*=-1
            IS.E=IS.E0-IS.h*IS.M

        end
    end
end

function Metropolis!(IS::Ising;steps::Int = 1,save_interval::Int=0,conf_save_interval::Int=0)
    @assert steps ≥ 0
    @assert save_interval ≥ 0

    ostart = IS.steps == 0 ? 0 : IS.steps+1
    nspins=length(IS.conf)
    saveidx::Int=1
    if save_interval >0
        nsave = length(ostart:save_interval:IS.steps+steps)
        M = zeros(Float64, nsave)
        E = zeros(Float64, nsave)
        if ostart==IS.steps
            M[1]=IS.M/nspins
            E[1]=IS.E/nspins
            saveidx=2
        end
    end

    for _ in 1:steps
        Metropolis_sweep!(IS)
        IS.steps += 1
        if save_interval>0 && IS.steps % save_interval == 0
            M[saveidx]=IS.M/nspins
            E[saveidx]=IS.E/nspins
            saveidx += 1
        end
        # if t ∈ 1:save_interval:steps
        #     selectdim(σ_t, ndims(σ) + 1, cld(t, save_interval)) .= σ
        # end
    end

    return save_interval>0 ? (E, M) : nothing
end
