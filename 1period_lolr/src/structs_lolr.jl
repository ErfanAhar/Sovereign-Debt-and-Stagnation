using LinearAlgebra
using QuantEcon

Base.@kwdef struct Model
    # Grid sizes
    Nb::Int = 1000          # private debt grid size
    Nl::Int = 100          # LOLR debt grid size
    Ng::Int = 2            # growth-state grid size
    Ne::Int = 10           # shock grid size

    # Grid bounds
    b_min::Float64 = -0.05
    b_max::Float64 = 1.0
    l_min::Float64 = 0.0
    l_max::Float64 = 0.0

    # Grids
    b::Vector{Float64} = collect(range(b_min, b_max, length = Nb)) # private debt grid
    l::Vector{Float64} = collect(range(l_min, l_max, length = Nl)) # LOLR debt grid
    g::Vector{Float64} = [0.96, 1.04]                              # growth states
    eps::Vector{Float64} = collect(range(-0.1, 0.1, length = Ne))  # shock grid

    # Transitions
    P_g::Matrix{Float64} = [0.60 0.40; 0.25 0.75]    # growth transition
    pi_eps::Vector{Float64} = fill(1.0 / Ne, Ne)     # iid shock probabilities

    # Preferences and technology
    beta::Float64 = 0.55                 # discount factor
    gamma::Float64 = 2.0                 # CRRA risk aversion
    sigma_eps::Float64 = 0.023           # shock scale in output
    n_std::Float64 = 2.5                 # Tauchen width
    rstar::Float64 = 0.01                # world interest rate
    theta::Float64 = 0.1                 # re-entry probability
    kappa::Float64 = 0.75                # debt haircut in default
    phi_g::Vector{Float64} = [0.985, 0.90] # output loss in default
    nbar_g::Vector{Float64} = [1.0, 1.0] # issuance cap

    # LOLR policy
    R_l_nd::Float64 = 1.05               # LOLR gross rate in repayment
    R_l_d::Float64 = 1.05                # LOLR gross rate in default
    Δ_nd::Float64 = 0.02                 # LOLR limit factor in repayment
    Δ_d::Float64 = 0.02                  # LOLR limit factor in default
    lbar_nd::Float64 = 0.02              # LOLR borrowing limit in repayment
    lbar_d::Float64 = 0.02               # LOLR borrowing limit in default

    # Schedule constraint
    pub::Float64 = 0.65                  # default-prob cutoff

    # Solver options
    max_iter::Int = 200                  # outer iterations
    max_iter_vd::Int = 200               # default value iterations
    max_iter_x::Int = 2000               # price iterations
    tol::Float64 = 1e-6                  # outer tolerance
    tol_vd::Float64 = 1e-6               # default tolerance
    tol_x::Float64 = 1e-6                # price tolerance
end

function init_model(m::Model = Model())

    b = m.b
    l = m.l
    g = m.g
    P_g = m.P_g

    Nb = length(b)
    Nl = length(l)
    Ng = length(g)
    Ne = m.Ne

    sigma_u = 1.0
    mc = QuantEcon.tauchen(Ne, 0.0, sigma_u, 0.0, m.n_std)
    eps = collect(mc.state_values)
    P_eps = mc.p
    P_eps_pow = P_eps^10000
    pi_eps = vec(P_eps_pow[1, :])
    pi_eps ./= sum(pi_eps)

    pi_g = fill(1.0 / Ng, Ng)
    for _ in 1:10_000
        pi_g = vec(pi_g' * P_g)
    end
    g_bar = dot(pi_g, g)

    lbar_nd_val = m.Δ_nd * g_bar
    lbar_d_val = m.Δ_d * g_bar
    l_max_val = max(m.l_max, lbar_nd_val, lbar_d_val)
    l = collect(range(m.l_min, l_max_val, length = Nl))

    fields = fieldnames(Model)
    base = NamedTuple{fields}(Tuple(getfield(m, f) for f in fields))
    return Model(; base..., Nb = Nb, Nl = Nl, Ng = Ng, Ne = Ne, b = b, l = l, g = g,
        eps = eps, P_g = P_g, pi_eps = pi_eps, l_max = l_max_val,
        lbar_nd = lbar_nd_val, lbar_d = lbar_d_val)
end

struct Solution
    # Value functions
    vnd::Array{Float64,4}       # value in repayment
    vd::Array{Float64,4}        # value in default

    # Default and re-entry
    d::BitArray{4}              # default indicator
    e::BitArray{4}              # re-entry indicator

    # Prices and schedules
    Q::Array{Float64,4}         # bond price
    X::Array{Float64,4}         # continuation price
    R::Array{Float64,3}         # gross rate schedule
    n::Array{Float64,3}         # issuance schedule

    # Default probabilities
    pdefault::Array{Float64,3}  # default probability

    # Policy indexes
    b_policy_idx::Array{Int,4}  # optimal b' index
    l_policy_idx::Array{Int,4}  # optimal l' index in repayment
    l_policy_idx_d::Array{Int,4} # optimal l' index in default

    # Diagnostics
    outer_iters::Int            # number of outer iterations
    vd_iters::Vector{Int}       # default-value iterations per outer loop
    x_iters::Vector{Int}        # price iterations per outer loop
    outer_errs::Vector{Float64} # vnd error per outer loop
    vd_errs::Vector{Float64}    # vd error per outer loop
    x_errs::Vector{Float64}     # price error per outer loop
    outer_times::Vector{Float64} # seconds per outer loop
end
