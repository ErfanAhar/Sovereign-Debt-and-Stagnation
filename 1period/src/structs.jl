using LinearAlgebra
using QuantEcon

Base.@kwdef struct Model
    # Grid sizes
    Nb::Int              # debt grid size
    Ng::Int               # growth-state grid size
    Ne::Int               # shock grid size

    # Grids
    b::Vector{Float64} = collect(range(-0.05, 0.55, length = Nb)) # debt grid
    g::Vector{Float64} = [0.96, 1.04]                             # growth states
    eps::Vector{Float64} = collect(range(-0.1, 0.1, length = Ne)) # shock grid

    # Transitions
    P_g::Matrix{Float64} = [0.60 0.40; 0.25 0.75]    # growth transition
    pi_eps::Vector{Float64} = fill(1.0 / Ne, Ne)     # iid shock probabilities

    # Preferences and technology
    beta::Float64 = 0.55                 # discount factor
    gamma::Float64 = 2.0                 # CRRA risk aversion
    sigma_eps::Float64 = 0.023           # shock scale in output
    rstar::Float64 = 0.01                # world interest rate
    theta::Float64 = 0.1                 # re-entry probability
    kappa::Float64 = 0.75                # debt haircut in default
    phi_g::Vector{Float64} = [0.985, 0.90] # output loss in default
    nbar_g::Vector{Float64} = [1.0, 1.0] # issuance cap

    # IMF policy
    R_l::Float64 = 1.02                            # IMF gross rate
    l_policy::Array{Float64,2} = fill(0.0, Ng, Ne)  # IMF loan rule

    # Schedule constraint
    pub::Float64 = 0.65                # default-prob cutoff

    # Solver options
    max_iter::Int = 200               # outer iterations
    max_iter_vd::Int = 200            # default value iterations
    max_iter_x::Int = 2000             # price iterations
    tol::Float64 = 1e-6               # outer tolerance
    tol_vd::Float64 = 1e-6            # default tolerance
    tol_x::Float64 = 1e-6             # price tolerance
end

function init_model(; Nb::Int = 3000, Ne::Int = 10, b_min::Float64 = -0.05, b_max::Float64 = 0.85,
    g::Vector{Float64} = [0.96, 1.04], P_g::Matrix{Float64} = [0.60 0.40; 0.25 0.75],
    sigma_eps::Float64 = 0.023, n_std::Float64 = 2.5,
    L0::Float64 = 0.0, l_policy::Union{Nothing, Array{Float64,2}} = nothing, kwargs...)

    Ng = length(g)
    b = collect(range(b_min, b_max, length = Nb))
    sigma_u = 1.0
    mc = QuantEcon.tauchen(Ne, 0.0, sigma_u, 0.0, n_std)
    eps = collect(mc.state_values)
    P_eps = mc.p
    P_eps_pow = P_eps^10000
    pi_eps = vec(P_eps_pow[1, :])
    pi_eps ./= sum(pi_eps)
    if l_policy === nothing
        l_policy = fill(L0, Ng, Ne)
    end

    return Model(; Nb = Nb, Ng = Ng, Ne = Ne, b = b, g = g, eps = eps, P_g = P_g,
        pi_eps = pi_eps, sigma_eps = sigma_eps, l_policy = l_policy, kwargs...)
end

struct Solution
    # Value functions
    vnd::Array{Float64,3}       # value in repayment
    vd::Array{Float64,3}        # value in default

    # Default and re-entry
    d::BitArray{3}              # default indicator
    e::BitArray{3}              # re-entry indicator

    # Prices and schedules
    Q::Array{Float64,3}         # bond price
    X::Array{Float64,3}         # continuation price
    R::Array{Float64,2}         # gross rate schedule
    n::Array{Float64,2}         # issuance schedule

    # Default probabilities
    pdefault::Array{Float64,2}  # default probability

    # Policy indexes
    b_policy_idx::Array{Int,3}  # optimal b' index

    # Diagnostics
    outer_iters::Int            # number of outer iterations
    vd_iters::Vector{Int}       # default-value iterations per outer loop
    x_iters::Vector{Int}        # price iterations per outer loop
    outer_errs::Vector{Float64} # vnd error per outer loop
    vd_errs::Vector{Float64}    # vd error per outer loop
    x_errs::Vector{Float64}     # price error per outer loop
    outer_times::Vector{Float64} # seconds per outer loop
end
