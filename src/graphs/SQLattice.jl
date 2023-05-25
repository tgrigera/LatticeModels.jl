# SQLattice.jl -- square lattice with different boundary conditions
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

export Site, connectivity, foreach_neighbour, foreach_neighbour!, foreach_bond
export SQLattice_open, SQLattice_periodic, SQLattice_site, random_site, random_site!
import Base.similar

#consider: abstract type Graph{T} <: AbstractArray{T} end

abstract type SQLattice{T} <: DenseArray{T,2} end

connectivity(l::SQLattice{T}) where T = 4

Base.size(lat::SQLattice{T}) where T = size(lat.nodes)

Base.getindex(lat::SQLattice{T},i::Int) where T = getindex(lat.nodes,i)

Base.getindex(lat::SQLattice{T},I::CartesianIndex{2}) where T = getindex(lat.nodes,I)

Base.setindex!(lat::SQLattice{T},v,i::Int) where T = setindex!(lat.nodes,v,i)

Base.setindex!(lat::SQLattice{T},v,I::CartesianIndex{2}) where T = setindex!(lat.nodes,v,I)

Base.IndexStyle(::Type{<:SQLattice}) = IndexLinear()

mutable struct SQLattice_site{Lat} <: Site where Lat <: SQLattice{T} where T
    nodes::Lat
    I::CartesianIndex{2}
end

SQLattice_site(lat::Lat, i::Int, j::Int) where Lat <: SQLattice{T} where T =
    SQLattice_site{Lat}(lat,CartesianIndex(i,j))


Base.getindex(lat::SQLattice{T},I::SQLattice_site{Lat}) where Lat <: SQLattice{T} where T = getindex(lat.nodes,I.I)

Base.setindex!(lat::SQLattice{T},v,I::SQLattice_site{Lat}) where Lat <: SQLattice{T} where T = setindex!(lat.nodes,v,I.I)

random_site(lat::Lat) where Lat <: SQLattice{T} where T = SQLattice_site(lat,rand(1:size(lat,1)),rand(1:size(lat,2)))

random_site(rng,lat::Lat) where Lat <: SQLattice{T} where T = SQLattice_site(lat,rand(rng,1:size(lat,1)),rand(rng,1:size(lat,2)))

function random_site!(site::SQLattice_site{Lat}) where Lat <: SQLattice{T} where T
    site.I=CartesianIndex(rand(1:size(site.nodes,1)),rand(1:size(site.nodes,1)))
end

function random_site!(rng,site::SQLattice_site{Lat}) where Lat <: SQLattice{T} where T
#    site.I=rand(rng,CartesianIndices(site.nodes))  # slower?
    site.I=CartesianIndex(rand(rng,1:size(site.nodes,1)),rand(rng,1:size(site.nodes,1)))
end

###############################################################################
#
# Open boundary conditions

struct SQLattice_open{T} <: SQLattice{T}
    nodes::Matrix{T}
end

SQLattice_open{T}(Lx,Ly) where T = SQLattice_open{T}(Matrix{T}(undef,Lx,Ly))

Base.similar(lat::SQLattice_open{T} where T,type::Type{N} where N) =
        SQLattice_open{type}(size(lat)...)

function foreach_neighbour(f,site::SQLattice_site{SQLattice_open{T}}) where T
    i,j=site.I[1],site.I[2]
    if i>1 f(site.nodes[i-1,j]) end
    if i<size(site.nodes,1) f(site.nodes[i+1,j]) end
    if j>1 f(site.nodes[i,j-1]) end
    if j<size(site.nodes,2) f(site.nodes[i,j+1]) end
end

function foreach_neighbour!(f,site::SQLattice_site{SQLattice_open{T}}) where T
    i,j=site.I[1],site.I[2]
    li=LinearIndices(site.nodes)
    if i>1 f(Ref(site.nodes,li[i-1,j])) end
    if i<size(site.nodes,1) f(Ref(site.nodes,li[i+1,j])) end
    if j>1 f(Ref(site.nodes,li[i,j-1])) end
    if j<size(site.nodes,2) f(Ref(site.nodes,li[i,j+1])) end
end

#
# Note that the is automatically column-major if using CartesianIndices
#
# This seems to be slower by a factor about 2
#
# function foreach_bond(f,lat::SQLattice_open{T}) where T
#     R=CartesianIndices(lat.nodes)
#     for I in first(R):CartesianIndex(last(R)[1]-1,last(R)[2])
#         f(lat[I],lat[I+CartesianIndex(1,0)])
#     end
#     for I in first(R):CartesianIndex(last(R)[1],last(R)[2]-1)
#         f(lat[I],lat[I+CartesianIndex(0,1)])
#     end
# end
function foreach_bond(f,lat::SQLattice_open{T}) where T
    L,M=size(lat.nodes)
    for j=1:M,i=1:L-1
        f(lat.nodes[i,j],lat.nodes[i+1,j])
    end
    for j=1:M-1,i=1:L
        f(lat.nodes[i,j],lat.nodes[i,j+1])
    end
end

#
# Distance binning
#

import BioStatPhys
import LinearAlgebra

function BioStatPhys.distance_binning(lat::SQLattice_open{T},Δr=1.) where T
    TII = Tuple{Int,Int}
    rmax = sqrt( (size(lat,1)-1)^2 + (size(lat,2)-1)^2 )
    binning = BioStatPhys.DistanceBinning(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,
                                          init=(x,n)->[TII[] for _=1:n])
    Indices=CartesianIndices(lat)
    LIndices=LinearIndices(Indices)
    for I ∈ Indices
        i = LIndices[I]
        for J ∈ Indices
            j = LIndices[J]
            push!(binning[LinearAlgebra.norm( Tuple(I-J) )],(i,j))
        end
    end
    return binning
end


###############################################################################
#
# Periodic boundary conditions

struct SQLattice_periodic{T} <: SQLattice{T}
    nodes::Matrix{T}
    N::Vector{Int}
    S::Vector{Int}
    E::Vector{Int}
    W::Vector{Int}
end

function SQLattice_periodic{T}(Lx,Ly) where T
    lat=SQLattice_periodic{T}(Matrix{T}(undef,Lx,Ly),Vector{Int}(undef,Ly),Vector{Int}(undef,Ly),
                              Vector{Int}(undef,Lx),Vector{Int}(undef,Lx))
    for i=2:Ly-1
        lat.N[i]=i-1
        lat.S[i]=i+1
    end
    lat.N[1]=Ly
    lat.N[Ly]=Ly-1
    lat.S[1]=2
    lat.S[Ly]=1
    for i=2:Lx-1
        lat.E[i]=i-1
        lat.W[i]=i+1
    end
    lat.E[1]=Lx
    lat.E[Lx]=Lx-1
    lat.W[1]=2
    lat.W[Lx]=1
    return lat
end

Base.similar(lat::SQLattice_periodic{T} where T,type::Type{N} where N) =
        SQLattice_periodic{type}(size(lat)...)

function foreach_neighbour(f,site::SQLattice_site{SQLattice_periodic{T}}) where T
    f(site.nodes[site.nodes.E[site.I[1]],site.I[2]])
    f(site.nodes[site.nodes.W[site.I[1]],site.I[2]])
    f(site.nodes[site.I[1],site.nodes.N[site.I[2]]])
    f(site.nodes[site.I[1],site.nodes.S[site.I[2]]])
end

neighbours(site::SQLattice_site{SQLattice_periodic{T}}) where T =
    [ site.nodes[site.nodes.E[site.I[1]],site.I[2]], 
      site.nodes[site.nodes.W[site.I[1]],site.I[2]], 
      site.nodes[site.I[1],site.nodes.N[site.I[2]]],
      site.nodes[site.I[1],site.nodes.S[site.I[2]]] ]

neighbour_indices(site::SQLattice_site{SQLattice_periodic{T}}) where T =
    CartesianIndex.([ (site.nodes.E[site.I[1]],site.I[2]), 
      (site.nodes.W[site.I[1]],site.I[2]), 
      (site.I[1],site.nodes.N[site.I[2]]),
      (site.I[1],site.nodes.S[site.I[2]]) ])

function foreach_neighbour!(f,site::SQLattice_site{SQLattice_periodic{T}}) where T
    li=LinearIndices(site.nodes)
    i,j=site.I[1],site.I[2]
    f(Ref(site.nodes,li[site.nodes.E[i],j]))
    f(Ref(site.nodes,li[site.nodes.W[i],j]))
    f(Ref(site.nodes,li[i,site.nodes.N[j]]))
    f(Ref(site.nodes,li[i,site.nodes.S[j]]))
end

function foreach_bond(f,lat::SQLattice_periodic{T}) where T
    L,M=size(lat.nodes)
    @inbounds for j=1:M,i=1:L-1
        f(lat.nodes[i,j],lat.nodes[i+1,j])
    end
    @inbounds for j=1:M
        f(lat.nodes[L,j],lat.nodes[1,j])
    end
    @inbounds for j=1:M-1,i=1:L
        f(lat.nodes[i,j],lat.nodes[i,j+1])
    end
    @inbounds for i=1:L
        f(lat.nodes[i,M],lat.nodes[i,1])
    end
end

@inline function ddiff(a::Int,b::Int,box_length::Int)
  temp=a-b
  return temp-box_length*round(Int,temp/box_length);
end

@inline function distancesq(lat::SQLattice_periodic{T},I,J) where T
    dx = ddiff(I[1],J[1],size(lat,1))
    dy = ddiff(I[2],J[2],size(lat,2))
   return dx*dx + dy*dy
end


function BioStatPhys.distance_binning(lat::SQLattice_periodic{T},Δr=1.) where T
    TII = Tuple{Int,Int}
    rmax = sqrt( (size(lat,1)-1)^2 + (size(lat,2)-1)^2 ) / 2
    binning=BioStatPhys.DistanceBinning(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,
                                       init=(x,n)->[TII[] for _=1:n])
    Indices=CartesianIndices(lat)
    LIndices=LinearIndices(Indices)
    for I ∈ Indices
        i=LIndices[I]
        for J ∈ Indices
            j=LIndices[J]
            push!(binning[ sqrt(distancesq(lat,I,J)) ],(i,j))
        end
    end
    return binning
end
