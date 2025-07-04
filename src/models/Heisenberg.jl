using Random
using StaticArrays
using LinearAlgebra
using JLD2
import Printf

export QHeisenberg, Metropolis_step!, update_observables!

"""
    QHeisenberg

Struct to hold a q-Heisenberg configuration plus parameters.  Energy
is

    J \\sum_{ij} n_{ij} (S_i-S_j)^2p - h \\cdot \\sum_i S_i

J is fixed to 1.

Alter only through `set_...` functions to ensure consistency.
"""
mutable struct QHeisenberg{q,N,Graph,RNG_T<:AbstractRNG}
    h::SVector{N,Float64}        # Magnetic field
    T::Float64                   # Temperature
    steps::Int64

    β::Float64                   # No J, so β is actually βJ
    M::MVector{N,Float64}        # Magnetisation
    m::Float64                   # Scalar magnetisaiton per spin
    E::Float64                   # Energy

    σ::Graph                     # spins
    RNG::RNG_T
end

"""
    QHeisenberg(graph_type,L...;seed=33,T=1.,h=0.,ordered=false)

Return a `QHeisenberg` object with underlying graph `graph_type`,
temperature `T` and field `h` (in the x direction).  A local random number generator with
seed `seed` will be created.  Spins will be will be set to a random
valuse or, if `ordered==true`, to 1.
"""
function QHeisenberg(graph_type,L...;q,N,seed=33,T=1.,h=0.,ordered=false)
    RNG = Random.Xoshiro(seed)
    graph = graph_type{MVector{N,Float64}}(L...)
    hv = [0. for _=1:N]
    hv[1] = h
    # note typeof(graph) below: this is required so that a concrete, rather than an abstract,
    # type, is passed as parameter to the QHeisenberg struct.
    qH = QHeisenberg{q,N,typeof(graph),Random.Xoshiro}(hv,0.,0,0.,hv,0.,0.,graph,RNG)
    for i ∈ eachindex(qH.σ)
        qH.σ[i] = MVector{N,Float64}([0. for _=1:N])
    end
    set_temperature!(qH,T)
    if ordered
        map(x->x[1]=1.,qH.σ)
    else
        for σ ∈ qH.σ
            rand_sphere!(RNG,σ)
        end
    end
    set_energy_mag!(qH)
    return qH
end

function set_temperature!(qH::QHeisenberg,T)
    qH.T = T
    qH.β = 1. /T
end    

function set_energy_mag!(qH::QHeisenberg)
    qH.E = energy_h0(qH)
    qH.M .= 0.
    for s ∈ qH.σ
        qH.M .+= s
    end
    qH.m = LinearAlgebra.norm(qH.M)/length(qH.σ)
    qH.E -= qH.h ⋅ qH.M
end

function energy_h0(qH::QHeisenberg{q}) where q
    E=Ref(0.)
    foreach_bond(qH.σ) do a,b
        E[] += (1-a⋅b)^q
    end
    return (2^q)*E[]
end

"""
    rand_sphere!(RNG,v::MVector{D,T}) where T<:Real

Store a random direction in D-dimensional Euclidean space (random
vector on the unit (D-1)-sphere) in `v`.  The algorithm for generic
`D` uses the method of forming a vector with D normally-distributed
(Gaussian with zero mean and unit variance) random numbers and then
normalizing.

There are specialized versions for D=2 and D=3.
"""
function rand_sphere!(RNG,v::MVector{D,T}) where D where T<:Real
    m = 0
    while true
        randn!(RNG,v)
        m = LinearAlgebra.norm(v)
        m==0. || break
    end
    v ./= m
end

# This method appears to be slower than the D-dimensional Gaussian method
# function rand_sphere!(RNG,v::MVector{2,T}) where T<:Real
#     θ = 2π * rand(RNG)
#     v[2],v[1] = sincos(θ)
# end


# This is from GSL, but it seems to run slower than the Gaussian one
# function rand_sphere!(RNG,r::MVector{2,T}) where T<:Real
#     # This method avoids trig, but it does take an average of 8/pi =
#     # 2.55 calls to the RNG, instead of one for the direct
#     # trigonometric method.
#     s = 0.
#     u = MVector(0.,0.)
#     while true
#         rand!(RNG,u)
#         u .= -1 .+ 2. * u
#         s = u[1] * u[1] + u[2] * u[2]
#         s > 1.0 || s == 0.0 || break
#     end

#     # This is the Von Neumann trick. See Knuth, v2, 3rd ed, p140
#     # (exercise 23).  Note, no sin, cos, or sqrt
#     r[1] = (u[1] * u[1] - u[2] * u[2]) / s
#     r[2] = 2 * u[1] * u[2] / s
# end

"""
    rand_sphere!(RNG,v::MVector{3,T}) where T<:Real

Store a random direction in 3-d (random vector on the unit 2-sphere)
in `v`.  This function uses a variant of the algorithm for computing a
random point on the unit sphere; the algorithm is suggested in Knuth,
v2, 3rd ed, p136; and attributed to Robert E Knop, CACM, 13 (1970),
326.

This coud be faster than the general method with Gaussian
distributions, depending on several factors.  On my laptop currently
it's about the same, or marginally slower, than the generic routine.

Translated from the GNU Scientific Library routine in `C`
"""
function rand_sphere!(RNG,v::MVector{3,T}) where T<:Real
    # Begin with the polar method for getting x,y inside a unit circle
    s = 0.
    while true
        v[1] = -1 + 2 * rand(RNG)
        v[2] = -1 + 2 * rand(RNG)
        s = v[1]*v[1] + v[2]*v[2]
        s>1.0 || break
    end      
    v[3] = -1 + 2 * s    # z uniformly distributed from -1 to 1
    a = 2 * sqrt(1 - s) # factor to adjust x,y so that x^2+y^2 is equal to 1-z^2
    v[1] *= a
    v[2] *= a
end

###############################################################################
#
# Single spin-flip Metropolis

function Metropolis_step!(qH::QHeisenberg{q}) where q
   I = random_site(qH.σ)
   ΔE = Ref{Float64}(0.)
   σn = copy(qH.σ[I])
    for _ ∈ 1:length(qH.σ)
        random_site!(qH.RNG,I)
        σo = qH.σ[I]
        ΔE[] = 0
        rand_sphere!(qH.RNG,σn)
        foreach_neighbour(I) do σ
            ΔE[] += (1 - σn⋅σ)^q - (1 - σo⋅σ)^q - qH.h⋅(σn .- σo)
        end
        ΔE[] *= 2^q
        if ΔE[]<=0. || rand(qH.RNG)<exp(-qH.β*ΔE[])
            qH.M .+= σn .- σo
            σo .= σn
            qH.E += ΔE[]
        end
    end
    qH.steps += 1
end

###############################################################################
##
## Log

function log_start(qH::QHeisenberg{q,N}) where q where N
    Ns = length(qH.σ)
    @info "Running $(q)-Heisenberg model with $N-dimensional spins\n"
    @info  "    Step          E/N      |m|\n";
    @info Printf.@sprintf(" Initial   %10.3e %10.3e\n",
	qH.E/Ns,qH.m)
end

function log(qH::QHeisenberg)
    Ns = length(qH.σ)
    @info Printf.@sprintf("%8ld   %10.3e %10.3e\n",
	qH.steps,qH.E/Ns,qH.m)
end

###############################################################################
##
## Observables

function update_observables!(qH::QHeisenberg)
    set_energy_mag!(qH)
end

###############################################################################
##
## Observers

struct StandardObserver{ModelT,IOT<:IO}# <: Observer{ModelT}
    filename::String
    io::IOT
    model::ModelT

    function StandardObserver(model,fname=nothing)
        if isnothing(fname) return new{typeof(model),typeof(stdout)}("",stdout,model) end
        return new{typeof(model),IOStream}(fname,open(fname,"w"),model)
    end
end

function observation_start(o::StandardObserver{QHeisenberg{q,N,graph,RNG},IOT}) where {q,N,graph,RNG,IOT}
    Printf.@printf(o.io,"#   Step    Energy    Magnetization Magnetization (vector)\n")
end

function observe(o::StandardObserver{QHeisenberg{q,N,graph,RNG},IOT}) where {q,N,graph,RNG,IOT}
    qH = o.model
    Ns = length(qH.σ)
    if N==2
        Printf.@printf(o.io,"%8ld   %10.3e %13.3e %10.3e %10.3e\n",
	    qH.steps,qH.E/Ns,qH.m,qH.M[1]/Ns,qH.M[2]/Ns)
    elseif N==3
        Printf.@printf(o.io,"%8ld   %10.3e %13.3e %10.3e %10.3e %10.3e\n",
	    qH.steps,qH.E/Ns,qH.m,qH.M[1]/Ns,qH.M[2]/Ns,qH.M[3]/Ns)
    else
        @error "N=$N not supported"
    end
end

function observation_stop(o::StandardObserver)
    if o.io!=stdout close(o.io) end
end

# struct MemoryObserver{N} <: Observer where N
#     M::Vector{MVector{N,Float64}}        # Magnetisation
#     m::Vector{Float64}                   # Scalar magnetisaiton per spin
#     E::Vector{Float64}                   # Energy
#     t::Vector{Int}
# end

# MemoryObserver(::QHeisenberg{q,N}) where q where N =
#     MemoryObserver{N}(MVector{N,Float64}[], Float64[], Float64[], Int[])

# function observe(qH::QHeisenberg,o::MemoryObserver)
#     push!(o.t,qH.steps)
#     push!(o.E,qH.E)
#     push!(o.M,qH.M)
#     push!(o.m,qH.m)
# end

mutable struct TrajectoryObserver{ModelT<:QHeisenberg,JIOT}
    qH          ::ModelT
    N            ::Int64
    fname        ::String
    jio          ::JIOT
    obs_interval ::Int64
    ssteps       ::Vector{Int}

    function TrajectoryObserver(model::QHeisenberg{q,N},fname,obs_interval) where {q,N}
        jio = jldopen(fname,"w")
        return new{typeof(model),typeof(jio)}(model,N,fname,jio,obs_interval,[])
    end
end

function observation_start(o::TrajectoryObserver)
    o.jio["N"] = o.N
    return nothing
end

function observe(o::TrajectoryObserver)
    push!(o.ssteps,o.qH.steps)
    o.jio["s$(o.qH.steps)/σ"] = o.qH.σ
    return nothing
end

function observation_stop(o::TrajectoryObserver)
    o.jio["ssteps"] = o.ssteps
    close(o.jio)
end

###############################################################################
##
## saving configuration and state

"""
    save_configuration(qH::QHeisenberg,filename)

Save lattice (spins) in binary format
"""
function save_configuration(qH::QHeisenberg{Q,N},filename) where {Q,N}
    jldsave(filename;qH.σ,N)
end

function load_configuration!(qH::QHeisenberg{Q,N},filename) where {Q,N}
    dic = load(filename)
    @assert dic["N"] == N
    qH.σ = dic["σ"]
    set_energy_mag!(qH)
end

function save_state(qH::QHeisenberg{Q,N},filename) where {Q,N}
     jldsave(filename; Q = Q, N = N, h = qH.h[1], T  = qH.T, steps = qH.steps,
         σ = qH.σ, RNG = qH.RNG, graph_type = Base.typename(typeof(qH.σ)).wrapper)
end

function load_state(::QHeisenberg,filename)
    dic = load(filename)
    L = [3 for _=1:dic["N"]]
    qH = QHeisenberg(dic["graph_type"],L...;q=dic["Q"],N=dic["N"],T=dic["T"],
         h=dic["h"], ordered=true)
    qH.steps = dic["steps"]
    qH.RNG = dic["RNG"]
    qH.σ = dic["σ"]
    set_energy_mag!(qH)
    return qH
end
