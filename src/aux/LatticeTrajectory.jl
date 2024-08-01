# LatticeTrajectory.jl -- functions to save lattice configurations
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

import Mmap

"""
    mutable struct LatticeTrajectory

This struct and the functions that operate with it help to create a
"trajectory" file, containing many snapshots of the lattice
configurations, that may be later accessed easily through memory
mapping.  If the configurations are `d`-dimensional, then using 'mmap'
below will return a `d+1`-dimensional array, where the last dimension
is a "time" index (index to different configurations).

Creating the object with the constructor
`LatticeTrajectory(::String,A)` will also open the file for writing.
"""
mutable struct LatticeTrajectory
    f
    ndims::Int
    nrec::Int
    size
end

"""
    ArrayFile(name::String,A)

Open and initialize file `name` (destroyed if it exists), and return
`LatticeTrajectory` object which can be used to write to it.

"""
function LatticeTrajectory(name::String,A)
    LT = LatticeTrajectory(open(name,"w+"),ndims(A),0,size(A))
    write(LT.f,LT.nrec)
    write(LT.f,LT.ndims)
    for i=1:LT.ndims
        write(LT.f,size(A,i))
    end
    return LT
end

import Base.write

function Base.write(f::LatticeTrajectory,A)
    f.nrec+=1
    write(f.f,A)
end

import Base.close

function Base.close(f::LatticeTrajectory)
    seek(f.f,0)
    write(f.f,f.nrec)
    close(f.f)
end

function mmap(name::String,t::T) where T<:Type
    f = open(name,"r")
    nrec=read(f,Int)
    ndims=read(f,Int)
    d=zeros(Int,ndims+1)
    for i=1:ndims
        d[i]=read(f,Int)
    end
    d[ndims+1]=nrec
    A=Mmap.mmap(f,Array{eltype(t),ndims+1},Tuple(d))
    return A,f
end
