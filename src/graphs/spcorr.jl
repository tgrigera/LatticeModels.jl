# spcorr.jl -- Space correlations for lattices
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


using AbstractFFTs
import Statistics

function space_correlation_Cr_Ck_iso(latset::AbstractVector;connected=false)

    (rx, ry, Crxy), (kx, ky, Ck) = space_correlation_Cr_Ck(latset,connected=connected)
    r,Cr = to_iso_2d(rx,ry,Crxy)
    k,Ck = to_iso_2d(kx,ky,Ck)

    return r,Cr,k,Ck
    
end

function space_correlation_Cr_Ck_iso(lat::Matrix{T};connected=false) where T<:Number

    (rx, ry, Crxy) , (kx, ky, Ck) = space_correlation_Cr_Ck(lat,connected=connected)
    r,Cr = to_iso_2d(rx,ry,Crxy)
    k,Ck = to_iso_2d(kx,ky,Ck)

    return r,Cr,k,Ck
    
end

function space_correlation_Cr_Ck_iso(lat::AbstractArray{T,3}) where T

    (rx, ry, rz, Cr) , (kx, ky, kz, Ck) = space_correlation_Cr_Ck(lat)
    r,Cr = to_iso_3d(rx,ry,rz,Cr)
    k,Ck = to_iso_3d(kx,ky,kz,Ck)

    return r,Cr,k,Ck
    
end

function space_correlation_Cr_Ck(latset::AbstractVector;connected=false)

    CC = space_correlation_Cr_Ck(latset[1],dry_run=true)
    Cr = last(first(CC))
    rr = CC[1][1:length(CC[1])-1]
    Ck = last(last(CC))
    rk = CC[2][1:length(CC[1])-1]

    for lat in latset
        CC1 = space_correlation_Cr_Ck(lat,connected=connected)
        Cr .+= last(first(CC1))
        Ck .+= last(last(CC1))
    end
    return (rr...,Cr/size(latset,1)),(rk...,Ck/size(latset,1))

end

function space_correlation_Cr_Ck(lat::AbstractMatrix{T};dry_run=false,connected=false) where T<:Number
    if connected
        lat2 = zeros(Float64,size(lat)...)
        mean = Statistics.mean(lat)
        lat2 .= lat .- mean
    else
        lat2 = lat
    end
    kx,ky,Ck = space_correlation_Ck_fft(lat2,dry_run)
    Crxy = dry_run ? zeros(Float64,size(lat)) : real.(ifft(Ck))
    rx = size(lat,1) * fftfreq(size(lat,1))
    ry = size(lat,2) * fftfreq(size(lat,2))
    return (rx, ry, Crxy) , (kx, ky, Ck)
end

function space_correlation_Cr_Ck(lat::AbstractArray{T,3},dry_run=false) where T<:Number
    kx,ky,kz,Ck = space_correlation_Ck_fft(lat,dry_run)
    Crxy = dry_run ? zeros(Float64,size(lat)) : real.(ifft(Ck))
    rx = size(lat,1) * fftfreq(size(lat,1))
    ry = size(lat,2) * fftfreq(size(lat,2))
    rz = size(lat,3) * fftfreq(size(lat,3))
    return (rx, ry, rz, Crxy) , (kx, ky, kz, Ck)
end

function to_iso_2d(rx,ry,C)
    rmax = max(maximum(rx),maximum(ry))
    Δ = min(rx[2],ry[2])
    Np = BioStatPhys.BinnedVector{Int}(Δ=Δ, min=Δ/2,max=rmax,round_max=RoundUp,
                                       init=zeros)
    Ciso = BioStatPhys.BinnedVector{Float64}(Δ=Δ, min=Δ/2,max=rmax,round_max=RoundUp,
                                              init=zeros)

    for (ik,kx) ∈ enumerate(rx), (jk,ky) ∈ enumerate(ry)
        k = sqrt(kx^2 + ky^2)
        Np[k] +=1
        Ciso[k] += C[ik,jk]
    end

    C0 = C[1,1]
    Ciso[0] -= C0
    Np[0] -= 1

    return vcat(0,collect(range(Ciso))),vcat(C0,Ciso./Np)
end

function to_iso_3d(rx,ry,rz,C)
    rmax = max(maximum(rx),maximum(ry),maximum(rz))
    Δ = min(rx[2],ry[2],rz[2])
    Np = BioStatPhys.BinnedVector{Int}(Δ=Δ,min=Δ/2,max=rmax,round_max=RoundUp,
                                       init=zeros)
    Ciso = BioStatPhys.BinnedVector{Float64}(Δ=Δ,min=Δ/2,max=rmax,round_max=RoundUp,
        init=zeros)

    for I ∈ CartesianIndices(C)
        k = sqrt(rx[I[1]]^2 + ry[I[2]]^2 + rz[I[3]]^2)
        Np[k] +=1
        Ciso[k] += C[I]
    end
    C0 = C[1,1]
    Ciso[0] -= C0
    Np[0] -= 1
                                              
    return vcat(0,collect(range(Ciso))), vcat(C0,Ciso./Np)
end

function space_correlation_Ck_fft(lat::AbstractMatrix{T},dry_run=false) where T<:Number
    Lx = size(lat,1)
    Ly = size(lat,2)
    N = Lx*Ly
    kxr = 2π * fftfreq(Lx)
    kyr = 2π * fftfreq(Ly)
    if dry_run return kxr,kyr,zeros(Float64,size(lat)) end

    Ck = abs2.(fft(lat)) ./ N
    return kxr,kyr,Ck
end

function space_correlation_Ck_fft(lat::AbstractArray{T,3},dry_run=false) where T<:Number
    Lx = size(lat,1)
    Ly = size(lat,2)
    Lz = size(lat,3)
    N = Lx*Ly*Lz
    kxr = 2π * fftfreq(Lx)
    kyr = 2π * fftfreq(Ly)
    kzr = 2π * fftfreq(Lz)


    Ck = abs2.(fft(lat)) ./ N
    return kxr,kyr,kzr,Ck
end

###############################################################################
#
# FFT routines for non-peridic lattices (2D)
#

function space_correlation_Cr_Ck_nonperiodic_iso(latset::AbstractVector;connected=false,Δ=nothing)
    (rx,ry,Cr),(kx,ky,Ck) = space_correlation_Cr_Ck_nonperiodic(latset,connected=connected,Δ=Δ)
    r,Cr = LatticeModels.to_iso_2d(rx,ry,Cr)
    k,Ck = LatticeModels.to_iso_2d(kx,ky,Ck)
    return r,Cr,k,Ck
end

function space_correlation_Cr_Ck_nonperiodic(latset::AbstractVector;connected=false,Δ=nothing)
    Cr = zero(latset[1])
    rx,ry = zeros(size(Cr,1)),zeros(size(Cr,2))
    Lkx,Lky = size(Cr).÷2 .+1
    Ck = zeros(Lkx,Lky)
    kx,ky = zeros(Lkx),zeros(Lky)

    for lat in latset
        (rx,ry,CC),(kx,ky,CCk) = space_correlation_Cr_Ck_nonperiodic(lat,connected=connected)
        Cr .+= CC
        Ck .+= CCk
    end
    if !isnothing(Δ)
        rx .*= Δ[1]
        ry .*= Δ[2]
        kx ./= Δ[1]
        ky ./= Δ[2]
    end
    return (rx,ry,Cr./size(latset,1)),(kx,ky,Ck./size(latset,1))
end

function space_correlation_Cr_Ck_nonperiodic(lat::AbstractMatrix{T};
                                             connected=false) where T<:Number
    Lx = size(lat,1)
    Ly = size(lat,2)
    lat2 = zeros(T,2Lx,2Ly)
    if connected
        mean = Statistics.mean(lat)
        lat2[1:Lx,1:Ly] .= lat .- mean
    else
        lat2[1:Lx,1:Ly] .= lat
    end
    (rx,ry,Cr),(kx,ky,Ck) = space_correlation_Cr_Ck(lat2)
    Cr = 4*Lx*Ly*Cr[1:Lx,1:Ly]
    for i=1:Lx Cr[i,:] ./= (Lx-i+1) end
    for j=1:Ly Cr[:,j] ./= (Ly-j+1) end

    # In Ck, half the frequencies must be discarded because the actual Δk is
    # twice that computed for the dataset padded with zeros.  After that,
    # we discard the negative frequencies, since Ck is even in k
    Ck = 4*Ck[1:2:Lx+1,1:2:Ly+1]
    kx = kx[1:2:Lx+1]
    kx[end] *= -1
    ky = ky[1:2:Ly+1]
    ky[end] *= -1

    return (rx[1:Lx],ry[1:Ly],Cr),(kx,ky,Ck)
end

###############################################################################
#
# Non-FFT routines, provided as debugging checks for the fft ones.  Not exported
#

function space_correlation_Ck_iso_direct(latset::AbstractVector)

    kxr,kyr,Ck = space_correlation_Ck_direct(latset)
    return to_iso_2d(kxr,kyr,Ck)

end

function space_correlation_Ck_direct(latset::AbstractVector)

    randC = space_correlation_Ck_direct(latset[1],dry_run=true)
    Ck = last(randC)
    rk = randC[1:length(randC)-1]

    for lat in latset
        rC1 = space_correlation_Ck_direct(lat)
        Ck .+= last(rC1)
    end
    return rk...,Ck/size(latset,1)

end

function space_correlation_Ck_direct(lat::AbstractMatrix{T};dry_run=false) where T<:Number

    Ck = zeros(Float64,size(lat))
    Lx = size(lat,1)
    Ly = size(lat,2)
    N = Lx*Ly
    kxr = 2π * fftfreq(Lx)
    kyr = 2π * fftfreq(Ly)
    if dry_run return kxr,kyr,Ck end

    for (ik,kx) ∈ enumerate(kxr), (jk,ky) ∈ enumerate(kyr)
        lk = Complex(0)
        for I ∈ CartesianIndices(lat)
            lk += exp( -im*(kx*(I[1]-1) + ky*(I[2]-1)) ) * lat[I]
        end
        Ck[ik,jk] = abs2(lk)/N
    end

    return kxr,kyr,Ck

end

"""
    DFT(v::Vector{T}) where T<:Number

Compute the 1-d discrete Fourier transform of vector `v`, which must
be indexed starting from 1.  The discrete Fourier transform is defined
as

`` F_{k+1} = \\sum_{j=0}^{N-1} \\exp(-2\\pi i j k) v_{j+1} ``

This is intended mainly for debugging purposes, since it should give
the same result as an FFT routine.
"""
function DFT(v::Vector{T}) where T<:Number
    N = size(v,1)
    Δ = 2π/N

    F = zeros(ComplexF64,N)
    for i ∈ 0:N-1
        for k ∈ 0:N-1
            F[i+1] += exp( -im*Δ*i*k ) * v[k+1]
        end
    end
    return F
end
