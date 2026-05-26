using LinearAlgebra
using QuantEcon

# One-period debt model with sunspots and ECB-style secondary-market purchases.
# Unlike the direct-LOLR model, the sovereign has only one debt state: private debt b.

Base.@kwdef mutable struct Model
    # State convention: g = 1 is low growth, s = 1 is the bad sunspot.
    Nb::Int = 1000
    Ng::Int = 2
    Ns::Int = 2
    Ne::Int = 17

    # State grids
    b_min::Float64 = -0.05
    b_max::Float64 = 1.00
    b::Vector{Float64} = collect(range(b_min, b_max, length = Nb))
    g::Vector{Float64} = [0.96, 1.04]
    eps::Vector{Float64} = collect(range(-0.1, 0.1, length = Ne))

    # Exogenous transitions
    P_g::Matrix{Float64} = [0.60 0.40; 0.25 0.75]
    P_s::Matrix{Float64} = [0.25 0.75; 0.25 0.75]
    pi_eps::Vector{Float64} = fill(1.0 / Ne, Ne)
    K::Matrix{Float64} = zeros(0, 0) # kron(P_s, P_g), filled by init_model

    # Preferences and default technology
    beta::Float64 = 0.725
    gamma::Float64 = 2.0
    sigma_eps::Float64 = 0.023
    n_std::Float64 = 2.5
    rstar::Float64 = 0.035
    theta::Float64 = 0.10
    kappa::Float64 = 0.75
    phi_g::Vector{Float64} = [0.93, 0.88]
    nbar_g::Vector{Float64} = [1.0, 1.0]

    # ECB facility policy from Section 2.1 of the note.
    # In repayment, low growth, and a bad sunspot, lenders may sell to the ECB
    # with probability intervention_prob at price 1 / (1 + rstar + ecb_spread).
    ecb_spread::Float64 = 0.04
    intervention_prob::Float64 = 0.50

    # Schedule-selection cutoff
    pub::Float64 = 0.65

    # Solver controls
    max_iter::Int = 200
    max_iter_vd::Int = 200
    max_iter_x::Int = 2000
    tol::Float64 = 1e-6
    tol_vd::Float64 = 1e-6
    tol_x::Float64 = 1e-6
end

struct Solution
    # Four-dimensional arrays use order (b, g, s, eps).
    vnd::Array{Float64,4}
    vd::Array{Float64,4}
    d::BitArray{4}
    e::BitArray{4}
    Q::Array{Float64,4}
    X::Array{Float64,4}

    # Schedule arrays use order (bprime, g, s).
    H::Array{Float64,3}        # payoff from holding the private bond
    L::Array{Float64,3}        # payoff from entering the ECB facility queue
    R::Array{Float64,3}        # gross rate offered to the sovereign
    n::Array{Float64,3}        # issuance implied by the debt choice
    pdefault::Array{Float64,3}
    schedule_mask::BitArray{3}
    ecb_used::BitArray{3}      # lenders strictly prefer entering the facility

    # Policy index for next-period debt.
    b_policy_idx::Array{Int,4}

    # Convergence diagnostics.
    outer_errs::Vector{Float64}
    vd_errs::Vector{Float64}
    x_errs::Vector{Float64}
end
