# Simple one-period model with private debt, LOLR debt, and sunspots.
# State: (b, l, g, s, eps)

function _stationary_mean_g(g::AbstractVector{<:Real}, P_g::AbstractMatrix{<:Real}; iters::Int = 10_000)
    pi_g = fill(1.0 / length(g), length(g))
    for _ in 1:iters
        pi_g = vec(pi_g' * P_g)
    end
    return sum(pi_g .* g)
end

Base.@kwdef mutable struct Model
    # Grid sizes
    Nb::Int = 1000
    Nl::Int = 1000
    Ng::Int = 2
    Ns::Int = 2
    Ne::Int = 17

    # Grid bounds
    b_min::Float64 = -0.05
    b_max::Float64 = 1.00
    l_min::Float64 = 0.00
    l_max::Float64 = 0.06

    # Grids
    b::Vector{Float64} = collect(range(b_min, b_max, length = Nb))
    l::Vector{Float64} = collect(range(l_min, l_max, length = Nl))
    g::Vector{Float64} = [0.96, 1.04]
    eps::Vector{Float64} = collect(range(-0.1, 0.1, length = Ne)) # shock grid

    # Transitions
    P_g::Matrix{Float64} = [0.60 0.40; 0.25 0.75]
    P_s::Matrix{Float64} = [0.25 0.75; 0.25 0.75]
    pi_eps::Vector{Float64} = fill(1.0 / Ne, Ne)
    K::Matrix{Float64} = zeros(0, 0)

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

    # LOLR borrowing terms
    R_l_nd::Float64 = 1.08
    R_l_d::Float64 = 1.08
    Δ_nd::Float64 = 0.06
    Δ_d::Float64 = 0.06
    lbar_nd::Float64 = Δ_nd * _stationary_mean_g(g, P_g)
    lbar_d::Float64 = Δ_d * _stationary_mean_g(g, P_g)

    # Schedule-selection cutoff
    pub::Float64 = 0.65

    # Solver controls
    max_iter::Int = 25
    max_iter_vd::Int = 40
    max_iter_x::Int = 60
    tol::Float64 = 1e-5
    tol_vd::Float64 = 1e-6
    tol_x::Float64 = 1e-6
end

struct Solution
    # Array order:
    # 5D arrays use (b, l, g, s, eps)
    # 4D schedule arrays use (b', l', g, s)

    # Value functions
    vnd::Array{Float64,5}
    vd::Array{Float64,5}

    # Default and re-entry
    d::BitArray{5}
    e::BitArray{5}

    # Prices and schedules
    Q::Array{Float64,5}
    X::Array{Float64,5}
    R::Array{Float64,4}
    n::Array{Float64,4}
    pdefault::Array{Float64,4}
    schedule_mask::BitArray{4}

    # Policy indexes
    b_policy_idx::Array{Int,5}
    l_policy_idx::Array{Int,5}
    l_policy_idx_d::Array{Int,5}

    # Diagnostics
    outer_errs::Vector{Float64}
    vd_errs::Vector{Float64}
    x_errs::Vector{Float64}
end
