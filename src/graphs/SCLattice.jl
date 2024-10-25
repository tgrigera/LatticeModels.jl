# SCLattice.jl -- simple cubic lattice with different boundary conditions
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

export SCLattice_periodic

abstract type SCLattice{T} <: DenseArray{T,3} end

connectivity(l::SCLattice{T}) where T = 6

Base.size(lat::SCLattice{T}) where T = size(lat.nodes)

Base.getindex(lat::SCLattice{T},i::Int) where T = getindex(lat.nodes,i)

Base.getindex(lat::SCLattice{T},I::CartesianIndex{2}) where T = getindex(lat.nodes,I)

Base.setindex!(lat::SCLattice{T},v,i::Int) where T = setindex!(lat.nodes,v,i)

Base.setindex!(lat::SCLattice{T},v,I::CartesianIndex{2}) where T = setindex!(lat.nodes,v,I)

Base.IndexStyle(::Type{<:SCLattice}) = IndexLinear()

mutable struct SCLattice_site{Lat} <: Site where Lat <: SCLattice{T} where T
    nodes::Lat
    I::CartesianIndex{3}
end

SCLattice_site(lat::Lat, i::Int, j::Int, k::Int) where Lat <: SCLattice{T} where T =
    SCLattice_site{Lat}(lat,CartesianIndex(i,j,k))

Base.getindex(lat::SCLattice{T},I::SCLattice_site{Lat}) where Lat <: SCLattice{T} where T = getindex(lat.nodes,I.I)

Base.setindex!(lat::SCLattice{T},v,I::SCLattice_site{Lat}) where Lat <: SCLattice{T} where T = setindex!(lat.nodes,v,I.I)

random_site(lat::Lat) where Lat <: SCLattice{T} where T = SCLattice_site(lat,rand(1:size(lat,1)),rand(1:size(lat,2)),rand(1:size(lat,3)))

random_site(rng,lat::Lat) where Lat <: SCLattice{T} where T = SCLattice_site(lat,rand(rng,1:size(lat,1)),rand(rng,1:size(lat,2)),rand(rng,1:size(lat,3)))

function random_site!(site::SCLattice_site{Lat}) where Lat <: SCLattice{T} where T
    site.I = CartesianIndex(rand(1:size(site.nodes,1)),rand(1:size(site.nodes,2)),
        rand(1:size(site.nodes,3)))
end

function random_site!(rng,site::SCLattice_site{Lat}) where Lat <: SCLattice{T} where T
    site.I = CartesianIndex(rand(rng,1:size(site.nodes,1)),rand(rng,1:size(site.nodes,2)),
        rand(rng,1:size(site.nodes,3)))
end

##############################################################################
#
# Periodic boundaries

struct SCLattice_periodic{T} <: SCLattice{T}
    nodes::Array{T,3}
    N::Vector{Int}
    S::Vector{Int}
    E::Vector{Int}
    W::Vector{Int}
    U::Vector{Int}
    D::Vector{Int}
end

function SCLattice_periodic{T}(Lx,Ly,Lz) where T
    lat = SCLattice_periodic{T}(Array{T,3}(undef,Lx,Ly,Lz),
          Vector{Int}(undef,Ly),Vector{Int}(undef,Ly),
          Vector{Int}(undef,Lx),Vector{Int}(undef,Lx),
          Vector{Int}(undef,Lz),Vector{Int}(undef,Lz))

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

    for i=2:Lz-1
        lat.U[i]=i-1
        lat.D[i]=i+1
    end
    lat.U[1]=Lz
    lat.U[Lz]=Lz-1
    lat.D[1]=2
    lat.D[Lz]=1

    return lat
end

Base.similar(lat::SCLattice_periodic{T} where T,type::Type{N} where N) =
        SCLattice_periodic{type}(size(lat)...)

function foreach_neighbour(f,site::SCLattice_site{SCLattice_periodic{T}}) where T
    f(site.nodes[site.nodes.E[site.I[1]],site.I[2],site.I[3]])
    f(site.nodes[site.nodes.W[site.I[1]],site.I[2],site.I[3]])
    f(site.nodes[site.I[1],site.nodes.N[site.I[2]],site.I[3]])
    f(site.nodes[site.I[1],site.nodes.S[site.I[2]],site.I[3]])
    f(site.nodes[site.I[1],site.I[2],site.nodes.U[site.I[3]]])
    f(site.nodes[site.I[1],site.I[2],site.nodes.D[site.I[3]]])
end

neighbours(site::SCLattice_site{SCLattice_periodic{T}}) where T =
    [ site.nodes[site.nodes.E[site.I[1]],site.I[2],site.I[3]],
     site.nodes[site.nodes.W[site.I[1]],site.I[2],site.I[3]],
     site.nodes[site.I[1],site.nodes.N[site.I[2]],site.I[3]],
     site.nodes[site.I[1],site.nodes.S[site.I[2]],site.I[3]],
     site.nodes[site.I[1],site.I[2],site.nodes.U[site.I[3]]],
     site.nodes[site.I[1],site.I[2],site.nodes.D[site.I[3]]] ]
        
neighbour_indices(site::SCLattice_site{SCLattice_periodic{T}}) where T =
    CartesianIndex.([ 
     (site.nodes.E[site.I[1]],site.I[2],site.I[3]),
     (site.nodes.W[site.I[1]],site.I[2],site.I[3]),
     (site.I[1],site.nodes.N[site.I[2]],site.I[3]),
     (site.I[1],site.nodes.S[site.I[2]],site.I[3]),
     (site.I[1],site.I[2],site.nodes.U[site.I[3]]),
     (site.I[1],site.I[2],site.nodes.D[site.I[3]]) ] )

function foreach_neighbour!(f,site::SCLattice_site{SCLattice_periodic{T}}) where T
    li = LinearIndices(site.nodes)
    i,j,k = site.I[1],site.I[2],site.I[3]
    f(Ref(site.nodes,li[site.nodes.E[i],j,k]))
    f(Ref(site.nodes,li[site.nodes.W[i],j,k]))
    f(Ref(site.nodes,li[i,site.nodes.N[j],k]))
    f(Ref(site.nodes,li[i,site.nodes.S[j],k]))
    f(Ref(site.nodes,li[i,j,site.nodes.U[k]]))
    f(Ref(site.nodes,li[i,j,site.nodes.D[k]]))
end

function foreach_bond(f,lat::SCLattice_periodic{T}) where T
    L,M,N=size(lat.nodes)
    @inbounds for k=1:N,j=1:M,i=1:L-1
        f(lat.nodes[i,j,k],lat.nodes[i+1,j,k])
    end
    @inbounds for k=1:N,j=1:M
        f(lat.nodes[L,j,k],lat.nodes[1,j,k])
    end
    @inbounds for k=1:N,j=1:M-1,i=1:L
        f(lat.nodes[i,j,k],lat.nodes[i,j+1,k])
    end
    @inbounds for k=1:N,i=1:L
        f(lat.nodes[i,M,k],lat.nodes[i,1,k])
    end
    @inbounds for k=1:N-1,j=1:M,i=1:L
        f(lat.nodes[i,j,k],lat.nodes[i,j,k+1])
    end
    @inbounds for j=1:M,i=1:L
        f(lat.nodes[i,j,N],lat.nodes[i,j,1])
    end
end
