
using Random
using Lattices

export Gaussian,set_temperature,enermag,sweep,MC

mutable struct Gaussian{Graph}
    J::Float64
    g::Float64
    T::Float64
    δ::Float64
    steps::Int64

    β::Float64
    βJ::Float64
    βg::Float64
    M::Float64        # Magnetisation
    E::Float64        # Energy

    conf::Graph
    RNG::AbstractRNG
end

function Gaussian(graph;seed=33,J=1.,T=1.,g=1.,δ=1.)
    RNG=Random.MersenneTwister(seed)
    G=Gaussian(J,g,T,δ,0,0.,0.,0.,0.,0.,graph,RNG)
    set_temperature(G,T)
    rand!(G.conf.nodes)
    enermag(G)
    return G
end

function set_temperature(G::Gaussian,T)
    G.T=T
    G.β=1. /T
    G.βJ=G.β*G.J
    G.βg=G.β*G.g
end

function enermag(G::Gaussian{T}) where T
    EJ=Ref(0.)
    foreach_bond(G.conf) do a,b
        EJ[]+=(a-b)^2
    end
    G.E=G.J*EJ[] + G.g * sum(x->x^2,G.conf.nodes)
    G.M=sum(G.conf.nodes)
end

function local_E(G::Gaussian,site,ϕ)
    S=Ref{Float64}(0.)
    foreach_neighbour(site) do x
        S[]+=(x-ϕ)^2
    end
    return G.J*S[] + G.g*ϕ^2
end

function sweep(G::Gaussian)
    site=random_site(G.RNG,G.conf)
    for _=1:length(G.conf)
        site=random_site(G.RNG,G.conf)
        E0=local_E(G,site,G.conf[site])
        ϕ=G.conf[site] + G.δ*2*(rand(G.RNG)-0.5)
        E1=local_E(G,site,ϕ)
        p=exp(G.β*(E0-E1))

        if p>1. || rand(G.RNG)<p
            G.E+=E1-E0
            G.M+=ϕ-G.conf[site]
            G.conf[site]=ϕ
        end
    end
end

function MC(G::Gaussian,steps,step0=0,Δ=1)
    st=Int64[]
    E=Float64[]
    M=Float64[]
    push!(st,step0)
    push!(E,G.E)
    push!(M,G.M)
    for s=1:steps
        sweep(G)
        if s%Δ==0
            enermag(G)
            push!(st,s+step0)
            push!(E,G.E)
            push!(M,G.M)
        end
    end
    enermag(G)
    return st,E,M
end
