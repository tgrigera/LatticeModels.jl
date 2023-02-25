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

abstract type SCLattice{T} <: DenseArray{T,3} end

connectivity(l::SCLattice{T}) where T = 6

Base.size(lat::SCLattice{T}) where T = size(lat.nodes)

struct SCLattice_periodic{T} <: SCLattice{T}
    nodes::Matrix{T}
end

function for_each_bond(f,lat::SCLattice_periodic{T}) where T
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
